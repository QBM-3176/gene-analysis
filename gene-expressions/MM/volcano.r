############################################################
# DESeq2 pipeline + Heatmap for GSE226353
# Input: GSE226353_ex_count.csv (raw counts)
# Output: Results CSVs + heatmap "GSE226353_heatmap_topXX.png"
############################################################

# Install packages if needed
# install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")
# install.packages("pheatmap")

library(DESeq2)
library(apeglm)
library(pheatmap)

# setwd("~/Striki/quant.2")

############################################################
# 1. Load count matrix
############################################################

# Path to your file (adjust if needed)
count_file <- "GSE226353_ex_count.csv"

# File format: semicolon-separated, first column index, 2nd Geneid, 3rd Length, then samples # nolint
raw <- read.csv(count_file, sep = ";", header = TRUE, stringsAsFactors = FALSE)

# Extract gene IDs and counts (drop index and Length columns)
counts <- raw[, -(1:3)]  # Samples start from column 4

# Verify dimensions and first columns (should match file header)
dim(counts)
colnames(counts)  # WT_1 WT_2 WT_3 CKD_veh_1 ...

############################################################
# 2. Create sample metadata (colData)
############################################################

# Exact column names from your file
sample_names <- colnames(counts)

# Conditions based on file header: WT (3), CKD_veh (3), CKD_RASi (3), CKD_RASi/SGLT2i (3), CKD_RASi/SGLT2i/MRA (3)
condition <- factor(c(
  rep("WT", 3),
  rep("CKD_veh", 3),
  rep("CKD_RASi", 3),
  rep("CKD_RASi_SGLT2i", 3),
  rep("CKD_RASi_SGLT2i_MRA", 3)
), levels = c("WT", "CKD_veh", "CKD_RASi", "CKD_RASi_SGLT2i", "CKD_RASi_SGLT2i_MRA"))

coldata <- data.frame(
  row.names = sample_names,
  condition = condition
)

table(coldata$condition)

############################################################
# 3. Construct DESeqDataSet
############################################################

dds_input <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),  # Ensure integer counts
  colData   = coldata,
  design    = ~ condition
)

# Prefilter low-count genes
keep <- rowSums(counts(dds_input)) >= 10
dds <- dds_input[keep, ]

dds_input

############################################################
# 4. Run DESeq2
############################################################

dds <- DESeq(dds_input)

############################################################
# 6. VST transform
############################################################

vsd <- vst(dds, blind = FALSE)

# Get VST expression matrix
vsd_mat <- assay(vsd)

############################################################
# 9. Volcano plots for each comparison
############################################################

# Define comparisons of interest
comparisons_1 <- list(
  c("WT", "CKD_veh")
)

comparisons_2 <- list(
  c("CKD_veh", "CKD_RASi"),
  c("CKD_veh", "CKD_RASi_SGLT2i"),
  c("CKD_veh", "CKD_RASi_SGLT2i_MRA")
)

# First level significance thresholds
pval_trhresh_1 <- 0.05
LFC_thresh_1 <- 1

# Second level significance thresholds
pval_trhresh_2 <- 0.01
LFC_thresh_2 <- 5

process_comparisons <- function(comparisons, dds) {
    for (comp in comparisons) {
    # Extract results
    res <- results(dds, contrast = c("condition", comp[1], comp[2]))
    res_df <- as.data.frame(res)

    # Remove NA values
    res_df <- res_df[!is.na(res_df$padj), ]

        # Add significance categories
    res_df$significance <- "Not significant"
    res_df$significance[res_df$padj < pval_trhresh_1 & res_df$log2FoldChange > LFC_thresh_1] <- "Upregulated"
    res_df$significance[res_df$padj < pval_trhresh_1 & res_df$log2FoldChange < -LFC_thresh_1] <- "Downregulated"

    # Color mapping
    volcano_colors <- c("Not significant" = "grey",
                        "Upregulated" = "red",
                        "Downregulated" = "blue")

    # Create volcano plot
    filename <- paste0("Volcano_", comp[1], "_vs_", comp[2], ".jpg")
    jpeg(filename, width = 2400, height = 2100, res = 300)
    par(mar = c(5, 5, 4, 2) + 0.1)

    plot(res_df$log2FoldChange, -log10(res_df$padj),
        col = volcano_colors[res_df$significance],
        pch = 19,
        cex = 0.6,
        xlab = "Log2 Fold Change",
        ylab = "-Log10 Adjusted P-value",
        main = paste0("Volcano Plot: ", comp[1], " vs ", comp[2]))

    # Add significance thresholds
    abline(h = -log10(pval_trhresh_1), col = "black", lty = 2, lwd = 1)
    abline(v = c(-LFC_thresh_1, LFC_thresh_1), col = "black", lty = 2, lwd = 1)

    # Add more significance thresholds
    abline(h = -log10(pval_trhresh_2), col = "cyan", lty = 1, lwd = 2)
    abline(v = c(-LFC_thresh_2, LFC_thresh_2), col = "cyan", lty = 1, lwd = 2)

    # Add legend
    legend("topright",
            legend = c(paste("Upregulated (", sum(res_df$significance == "Upregulated"), ")", sep = ""),
                    paste("Downregulated (", sum(res_df$significance == "Downregulated"), ")", sep = ""),
                    "Not significant"),
            col = c("red", "blue", "grey"),
            pch = 19,
            cex = 0.8)

    dev.off()
    cat(paste0("Volcano plot saved as ", filename, "\n"))

    # Extract most significant downregulated genes by padj and effect size
    sign_res_df <- res_df[res_df$padj < pval_trhresh_2 & abs(res_df$log2FoldChange) > LFC_thresh_2, ]
    sign_res_df$Gene <- rownames(sign_res_df)
    
    # Reorder columns to put Gene first
    sign_res_df <- sign_res_df[, c("Gene", setdiff(names(sign_res_df), "Gene"))]

    # Save to CSV
    filename <- paste0("MSDEG_", comp[1], "_vs_", comp[2], ".csv")
    write.csv(sign_res_df, filename, row.names = FALSE)

    cat(paste0("Most significant differentially expressed genes saved as ", filename, "\n"))
    cat("\n")
  }
}

process_comparisons(comparisons_1, dds)
process_comparisons(comparisons_2, dds)

# Generate volcano plot for each comparison
# for (comp in comparisons_1) {
#   # Extract results
#   res <- results(dds, contrast = c("condition", comp[1], comp[2]))
#   res_df <- as.data.frame(res)

#   # Remove NA values
#   res_df <- res_df[!is.na(res_df$padj), ]

#     # Add significance categories
#   res_df$significance <- "Not significant"
#   res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
#   res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

#   # Color mapping
#   volcano_colors <- c("Not significant" = "grey",
#                       "Upregulated" = "red",
#                       "Downregulated" = "blue")

#   # Create volcano plot
#   filename <- paste0("QC_volcano_", comp[1], "_vs_", comp[2], ".jpg")
#   jpeg(filename, width = 2400, height = 2100, res = 300)
#   par(mar = c(5, 5, 4, 2) + 0.1)

#   plot(res_df$log2FoldChange, -log10(res_df$padj),
#        col = volcano_colors[res_df$significance],
#        pch = 19,
#        cex = 0.6,
#        xlab = "Log2 Fold Change",
#        ylab = "-Log10 Adjusted P-value",
#        main = paste0("Volcano Plot: ", comp[1], " vs ", comp[2]))

#   # Add significance thresholds
#   abline(h = -log10(0.05), col = "black", lty = 2, lwd = 1)
#   abline(v = c(-1, 1), col = "black", lty = 2, lwd = 1)

#   # Add legend
#   legend("topright",
#          legend = c(paste("Upregulated (", sum(res_df$significance == "Upregulated"), ")", sep = ""),
#                    paste("Downregulated (", sum(res_df$significance == "Downregulated"), ")", sep = ""),
#                    "Not significant"),
#          col = c("red", "blue", "grey"),
#          pch = 19,
#          cex = 0.8)

#   dev.off()

#   cat(paste0("Volcano plot saved as ", filename, "\n"))

#   # Extract top 10 most significant downregulated genes by padj and effect size
#   # Sort by padj (most significant first), then by absolute log2FoldChange
#   res_df_sorted_downregulated <- res_df[order(res_df$padj, -abs(res_df$log2FoldChange)), ]

#   # Get top 10 downregulated genes
#   top10_downregulated_genes <- head(res_df_sorted_downregulated, 10)

#   # Add gene names (rownames) as a column
#   top10_downregulated_genes$Gene <- rownames(top10_downregulated_genes)

#   # Reorder columns to put Gene first
#   top10_downregulated_genes <- top10_downregulated_genes[, c("Gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "significance")]

#   # Save to CSV
#   top10_downregulated_filename <- paste0("Top10_downregulated_genes_", comp[1], "_vs_", comp[2], ".csv")
#   write.csv(top10_downregulated_genes, top10_downregulated_filename, row.names = FALSE)

#   cat(paste0("Top 10 downregulated genes saved as ", top10_downregulated_filename, "\n"))
#   cat("Top 10 most significantly differentially expressed downregulated genes:\n")
#   print(top10_downregulated_genes[, c("Gene", "log2FoldChange", "padj")])
#   cat("\n")



#   # Extract top 10 most significant upregulated genes by padj and effect size
#   # Sort by padj (most significant first), then by absolute log2FoldChange
#   res_df_sorted_upregulated <- res_df[order(res_df$padj, -abs(res_df$log2FoldChange), decreasing = TRUE), ]

#   # Get top 10 upregulated genes
#   top10_upregulated_genes <- head(res_df_sorted_upregulated, 10)

#   # Add gene names (rownames) as a column
#   top10_upregulated_genes$Gene <- rownames(top10_upregulated_genes)
#   # Reorder columns to put Gene first
#   # top10_upregulated_genes <- top10_upregulated_genes[, c("Gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "significance")]

#   # Save to CSV
#   top10_upregulated_filename <- paste0("Top10_upregulated_genes_", comp[1], "_vs_", comp[2], ".csv")
#   write.csv(top10_upregulated_genes, top10_upregulated_filename, row.names = FALSE)

#   cat(paste0("Top 10 genes saved as ", top10_upregulated_filename, "\n"))
#   cat("Top 10 most significantly differentially expressed upregulated genes:\n")
#   print(top10_upregulated_genes[, c("Gene", "log2FoldChange", "padj")])
#   cat("\n")

# }
