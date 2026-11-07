# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# Install packages if needed
# install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("clusterProfiler")
# BiocManager::install("enrichplot")


library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)

set.seed(42)

# Load the data
counts <- read.csv("GSE226353_ex_count.csv", sep = ";", row.names = 1)
rownames(counts) <- counts$Geneid
gene_info <- counts[, 1:2] # Geneid and Length
count_data <- counts[, 3:ncol(counts)] # Just the 15 count columns

# Create metadata
col_data <- data.frame(
  sample = colnames(count_data),
  condition = factor(c(rep("WT", 3), rep("CKD_veh", 3), rep("CKD_RASi", 3), 
                       rep("CKD_SGLT2", 3), rep("CKD_Triple", 3)))
)
rownames(col_data) <- col_data$sample

# Setup DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = round(count_data),
                              colData = col_data,
                              design = ~ condition)

# Run Differential Expression
dds <- DESeq(dds)

# Get results for a specific comparison (e.g., Disease vs healthy)
res <- results(dds, contrast=c("condition", "CKD_veh", "WT"))

# Create a ranked list based on log2FoldChange
# We remove NAs and sort from highest change to lowest
res_df <- as.data.frame(res)
res_df$gene <- rownames(counts) # Using Ensembl IDs
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

# Save GSEA results table
gsea_table <- gsea_results@result
gsea_table <- gsea_table[, -ncol(gsea_table), drop = FALSE]
write.csv(gsea_table, "PathwayGSEA.csv", row.names = FALSE)

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

# Enrichment plot for a specific pathway (e.g., the top one)
mountain_plot <- gseaplot2(gsea_results, geneSetID = 1, title = gsea_results$Description[1])
ggsave("PathwayGSEA_CKD_veh_vs_WT-mountain.jpg", mountain_plot, device = "jpg", units = "cm", width = 16, height = 10)

dotplot <- dotplot(gsea_results, showCategory = 10, split = ".sign") + facet_grid(.~.sign) +
  scale_y_discrete(labels = truncate_desc) +
  theme(axis.text.y = element_text(size = 10))
ggsave("PathwayGSEA_CKD_veh_vs_WT-dotplot.jpg", dotplot, device = "jpg", units = "cm", width = 16, height = 16)

# --- Bar plot of GSEA results ---


gsea_df <- gsea_results@result
# Select top 10 activated and top 10 suppressed pathways by NES
top_activated <- gsea_df[gsea_df$NES > 0, ]
top_activated <- top_activated[order(-top_activated$NES), ][1:min(10, nrow(top_activated)), ]
top_suppressed <- gsea_df[gsea_df$NES < 0, ]
top_suppressed <- top_suppressed[order(top_suppressed$NES), ][1:min(10, nrow(top_suppressed)), ]
top_pathways <- rbind(top_activated, top_suppressed)
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

ggsave("PathwayGSEA_CKD_veh_vs_WT-barplot.jpg", device = "jpg", units = "cm", width = 16, height = 16)

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
print(p_up)
ggsave("PathwayGSEA_CKD_veh_vs_WT-barplot_up.jpg", p_up, device = "jpg", units = "cm", width = 16, height = 16)

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
ggsave("PathwayGSEA_CKD_veh_vs_WT-barplot_down.jpg", p_down, device = "jpg", units = "cm", width = 16, height = 16)