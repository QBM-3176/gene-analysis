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
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)

set.seed(42)

# setwd("~/Striki/quant.2")

############################################################
# 1. Load count matrix
############################################################

# Path to your file (adjust if needed)
count_file <- "GSE226353_ex_count.csv"

# Load the data
counts <- read.csv("GSE226353_ex_count.csv", sep = ";", row.names = 1)
rownames(counts) <- counts$Geneid
gene_info <- counts[, 1:2] # Geneid and Length
count_data <- counts[, 3:ncol(counts)] # Just the 15 count columns

# Verify dimensions and first columns (should match file header)
dim(counts)
colnames(counts)  # WT_1 WT_2 WT_3 CKD_veh_1 ...

############################################################
# 2. Create sample metadata (colData)
############################################################

# Exact column names from your file
sample_names <- colnames(count_data)

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
  countData = round(as.matrix(count_data)),  # Ensure integer counts
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
# VST transform
############################################################

vsd <- vst(dds, blind = FALSE)

# Get VST expression matrix
vsd_mat <- assay(vsd)

# Define comparisons of interest
comparisons_1 <- list(
  c("CKD_veh", "WT")
)

comparisons_2 <- list(
  c("CKD_RASi", "CKD_veh"),
  c("CKD_RASi_SGLT2i", "CKD_veh"),
  c("CKD_RASi_SGLT2i_MRA", "CKD_veh")
)

# Truncate descriptions to max 40 characters without cutting words
truncate_desc <- function(x, max_chars = 40) {
  sapply(x, function(s) {
    if (nchar(s) <= max_chars) return(s)
    words <- strsplit(s, " ")[[1]]
    result <- ""
    for (w in words) {
      candidate <- if (result == "") w else paste(result, w)
      if (nchar(candidate) > max_chars) break
      result <- candidate
    }
    paste0(result, "...")
  }, USE.NAMES = FALSE)
}

process_comparisons <- function(comparisons, dds) {
    for (comp in comparisons) {
    # Get results for a specific comparison (e.g., Disease vs healthy)
    res <- results(dds, contrast=c("condition", comp[1], comp[2]))

    # Create a ranked list based on log2FoldChange
    # We remove NAs and sort from highest change to lowest
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df) # Using Ensembl IDs
    res_df <- res_df[!is.na(res_df$log2FoldChange), ]

    gene_list <- res_df$log2FoldChange
    names(gene_list) <- res_df$gene
    gene_list <- sort(gene_list, decreasing = TRUE)

    # Run GSEA for Biological Processes (BP)
    gsea_results <- gseGO(geneList = gene_list,
                        OrgDb = org.Mm.eg.db,
                        keyType = "ENSEMBL",
                        ont = "BP",
                        pvalueCutoff = 0.05,
                        verbose = TRUE)

    # Convert IDs to readable gene symbols for the final table
    gsea_results <- setReadable(gsea_results, OrgDb = org.Mm.eg.db)

    # print(head(gsea_results@result))

    # Enrichment plot for a specific pathway (e.g., the top one)
    title <- paste0(gsea_results$Description[1], "  [", comp[1], " vs ", comp[2], "]")
    p <- gseaplot2(gsea_results, geneSetID = 1, title = title, base_size = 16)
    plot_filename <- paste0("PathwayGSEA_", comp[1], "_vs_", comp[2], ".jpg")
    ggsave(plot_filename, plot = p, width = 10, height = 6, dpi = 300)
    cat(paste0("GSEA plot saved as ", plot_filename, "\n"))

    # Save GSEA results table
    filename <- paste0("PathwayGSEA_", comp[1], "_vs_", comp[2])
    gsea_table <- gsea_results@result
    gsea_table <- gsea_table[, -ncol(gsea_table), drop = FALSE]
    write.csv(gsea_table, filename, row.names = FALSE)

    cat(paste0("Most significant differentially expressed pathways saved as ", filename, ".csv", "\n"))
    cat("\n")

    # ========================================================================


    # Enrichment plot for a specific pathway (e.g., the top one)
    # mountain_plot <- gseaplot2(gsea_results, geneSetID = 1, title = gsea_results$Description[1])
    # ggsave("PathwayGSEA_CKD_veh_vs_WT-mountain.jpg", mountain_plot, device = "jpg", units = "cm", width = 16, height = 10)

    filename_dotplot <- paste0(filename, "-dotplot.jpg")
    dotplot <- dotplot(gsea_results, showCategory = 10, split = ".sign") + facet_grid(.~.sign) +
      scale_y_discrete(labels = truncate_desc) +
      theme(axis.text.y = element_text(size = 10))
    ggsave(filename_dotplot, dotplot, device = "jpg", units = "cm", width = 16, height = 16)

    # --- Bar plot of GSEA results ---


    gsea_df <- gsea_results@result

    # Select top 10 activated and top 10 suppressed pathways by NES
    top_up <- gsea_df[gsea_df$NES > 0, ]
    top_up <- top_up[order(-top_up$NES), ][1:min(10, nrow(top_up)), ]
    top_down <- gsea_df[gsea_df$NES < 0, ]
    top_down <- top_down[order(top_down$NES), ][1:min(10, nrow(top_down)), ]

    top_pathways <- rbind(top_up, top_down)
    top_pathways <- top_pathways[order(top_pathways$NES), ]
    # Truncate long pathway names to max 40 characters
    top_pathways$Description <- truncate_desc(top_pathways$Description, 40)
    top_pathways$Description <- factor(top_pathways$Description, levels = top_pathways$Description)

    ggplot(top_pathways, aes(x = NES, y = Description, fill = p.adjust)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient(low = "red", high = "blue", name = "P-adjust") +
      theme_minimal(base_size = 12) +
      labs(x = "Normalized Enrichment Score (NES)", y = NULL,
          title = "GSEA: Top Enriched Pathways",
          subtitle = "CKD_veh vs WT") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
      theme(axis.text.y = element_text(size = 10))

    filename_barplot <- paste0(filename, "-barplot.jpg")
    ggsave(filename_barplot, device = "jpg", units = "cm", width = 16, height = 16)

    # --- Separate bar plots for top 10 up- and down-regulated pathways ---

    # Top 10 upregulated pathways
    top_up <- gsea_df[gsea_df$NES > 0, ]
    top_up <- top_up[order(-top_up$NES), ][1:min(10, nrow(top_up)), ]
    top_up$Description <- truncate_desc(top_up$Description, 40)
    top_up$Description <- factor(top_up$Description, levels = rev(top_up$Description))

    p_up <- ggplot(top_up, aes(x = NES, y = Description, fill = p.adjust)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient(low = "red", high = "blue", name = "P-adjust") +
      theme_minimal(base_size = 12) +
      labs(x = "Normalized Enrichment Score (NES)", y = NULL,
          title = "Top 10 Upregulated Pathways",
          subtitle = "CKD_veh vs WT") +
      scale_y_discrete(labels = truncate_desc) +
      theme(axis.text.y = element_text(size = 10))
    filename_barplot_up <- paste0(filename, "-barplot_up.jpg")
    ggsave(filename_barplot_up, p_up, device = "jpg", units = "cm", width = 16, height = 16)

    # Top 10 downregulated pathways
    top_down <- gsea_df[gsea_df$NES < 0, ]
    top_down <- top_down[order(top_down$NES), ][1:min(10, nrow(top_down)), ]
    top_down$Description <- truncate_desc(top_down$Description, 40)
    top_down$Description <- factor(top_down$Description, levels = rev(top_down$Description))

    p_down <- ggplot(top_down, aes(x = NES, y = Description, fill = p.adjust)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient(low = "red", high = "blue", name = "P-adjust") +
      theme_minimal(base_size = 12) +
      labs(x = "Normalized Enrichment Score (NES)", y = NULL,
          title = "Top 10 Downregulated Pathways",
          subtitle = "CKD_veh vs WT") +
      theme(axis.text.y = element_text(size = 10))
    print(p_down)
    filename_barplot_down <- paste0(filename, "-barplot_down.jpg")
    ggsave(filename_barplot_down, p_down, device = "jpg", units = "cm", width = 16, height = 16)
    }
}

process_comparisons(comparisons_1, dds)
process_comparisons(comparisons_2, dds)
