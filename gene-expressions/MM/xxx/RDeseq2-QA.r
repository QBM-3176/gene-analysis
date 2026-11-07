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
# 5. QC plots
############################################################

plotDispEsts(dds)
plotMA(results(dds, contrast = c("condition", "CKD_RASi_SGLT2i_MRA", "WT")), ylim = c(-8, 8))

############################################################
# 6. VST transform
############################################################

vsd <- vst(dds, blind = FALSE)

############################################################
# 7. Boxplots by condition group
############################################################

# Get VST expression matrix
vsd_mat <- assay(vsd)

# Prepare data for boxplot grouped by condition
boxplot_data <- data.frame(
  expression = as.vector(vsd_mat),
  sample = rep(colnames(vsd_mat), each = nrow(vsd_mat)),
  condition = rep(coldata$condition, each = nrow(vsd_mat))
)

# Define colors matching your annotation colors
group_colors <- c(
  WT = "#00FFFF",
  CKD_veh = "#0000FF",
  CKD_RASi = "#FF00FF",
  CKD_RASi_SGLT2i = "#00FF00",
  CKD_RASi_SGLT2i_MRA = "#FF0000"
)

# Create boxplots
jpeg("QC_boxplots_by_condition.jpg", width = 3000, height = 1800, res = 300)
par(mar = c(8, 4, 4, 2) + 0.1)
boxplot(expression ~ condition, data = boxplot_data,
        main = "VST Expression Distribution by Condition",
        ylab = "VST Expression",
        xlab = "",
        las = 2,
        col = group_colors[levels(boxplot_data$condition)],
        outline = FALSE)
dev.off()

cat("Boxplot saved as QC_boxplots_by_condition.jpg\n")

############################################################
# 8. PCA plot with custom colors
############################################################

# Perform PCA on VST-transformed data (transpose so samples are rows)
pca_result <- prcomp(t(vsd_mat), scale. = FALSE)

# Extract PC1 and PC2 scores
pc_scores <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  sample = colnames(vsd_mat),
  condition = coldata$condition
)

# Calculate variance explained by each PC
var_explained <- round(100 * summary(pca_result)$importance[2, 1:2], 1)

# Create PCA plot
jpeg("QC_PCA_colored.jpg", width = 2400, height = 2100, res = 300)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(pc_scores$PC1, pc_scores$PC2,
     col = group_colors[as.character(pc_scores$condition)],
     pch = 19,
     cex = 2,
     xlab = paste0("PC1 (", var_explained[1], "% variance)"),
     ylab = paste0("PC2 (", var_explained[2], "% variance)"),
     main = "PCA of VST-transformed Expression Data")

# Add legend
legend("topright",
       legend = levels(pc_scores$condition),
       col = group_colors[levels(pc_scores$condition)],
       pch = 19,
       cex = 0.8,
       pt.cex = 1.5)

# Add sample labels
text(pc_scores$PC1, pc_scores$PC2, labels = pc_scores$sample,
     pos = 3, cex = 0.6)

dev.off()

cat("PCA plot saved as QC_PCA_colored.jpg\n")

############################################################
# 9. Volcano plots for each comparison
############################################################

# Define comparisons of interest
comparisons <- list(
  c("CKD_veh", "WT"),
  c("CKD_RASi", "WT"),
  c("CKD_RASi_SGLT2i", "WT"),
  c("CKD_RASi_SGLT2i_MRA", "WT")
)

# Generate volcano plot for each comparison
for (comp in comparisons) {
  # Extract results
  res <- results(dds, contrast = c("condition", comp[1], comp[2]))
  res_df <- as.data.frame(res)

  # Remove NA values
  res_df <- res_df[!is.na(res_df$padj), ]

  # Add significance categories
  res_df$significance <- "Not significant"
  res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
  res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

  # Color mapping
  volcano_colors <- c("Not significant" = "grey", 
                      "Upregulated" = "red", 
                      "Downregulated" = "blue")

  # Create volcano plot
  filename <- paste0("QC_volcano_", comp[1], "_vs_", comp[2], ".jpg")
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
  abline(h = -log10(0.05), col = "black", lty = 2, lwd = 1)
  abline(v = c(-1, 1), col = "black", lty = 2, lwd = 1)

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
}

############################################################
# 10. Heatmap: Top TOP_GENES_COUNT most variable genes
############################################################

# Select top variable genes,
all_genes_sorted <- names(sort(apply(assay(vsd), 1, var), decreasing = TRUE))
# genes_with_symbols <- all_genes_sorted[!grepl("^ENSMUSG", all_genes_sorted)]
# top_var_genes <- genes_with_symbols[1:TOP_GENES_COUNT]
top_var_genes <- all_genes_sorted[1:TOP_GENES_COUNT]

# Subset VST matrix
mat <- assay(vsd)[top_var_genes, ]

# Transpose matrix (genes become columns, samples become rows)
mat <- t(mat)

# Row/column annotation colors (matching uploaded image style)
annotation_row <- data.frame(
  Condition = coldata$condition,
  row.names = rownames(mat)
)
annotation_colors <- list(
  Condition = c(
    WT = "#00FFFF",
    CKD_veh = "#0000FF",
    CKD_RASi = "#FF00FF",
    CKD_RASi_SGLT2i =  "#00FF00",
    CKD_RASi_SGLT2i_MRA =  "#FF0000"
  )
)

# Generate heatmap (matches uploaded image: red-blue, clustered, color bar 1.5 to -1.5)
# pheatmap(mat,
#          scale = "column",  # Z-score per gene (column-centered)
#          clustering_distance_rows = "euclidean",
#          clustering_distance_cols = "euclidean",
#          clustering_method = "complete",
#          color = colorRampPalette(c("navy", "white", "red"))(TOP_GENES_COUNT),
#          breaks = seq(-1.5, 1.5, length.out = TOP_GENES_COUNT+1),  # Exact -1.5 to +1.5 range
#          annotation_row = annotation_row,
#          annotation_colors = annotation_colors,
#          show_rownames = TRUE,  # Show sample labels
#          show_colnames = TRUE,  # Show gene labels
#          fontsize_col = 8,
#          main = paste0("GSE226353: Top ",  TOP_GENES_COUNT, " Variable Genes"),
#          filename = paste0("GSE226353_heatmap_top", TOP_GENES_COUNT, ".png"),  # Saves PNG
#          width = 12, height = 4)

# cat(paste0("Heatmap saved as GSE226353_heatmap_top", TOP_GENES_COUNT, ".png\n"))
