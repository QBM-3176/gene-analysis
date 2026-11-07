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
# BiocManager::install("org.Mm.eg.db")
# install.packages("ggplot2")
# install.packages("gridExtra")
# install.packages("reshape2")

library(DESeq2)
library(apeglm)
library(pheatmap)
library(org.Mm.eg.db)
library(ggplot2)
library(gridExtra)
library(reshape2)

############################################################
# 1. Load count matrix
############################################################

# Path to your file (adjust if needed)
count_file <- "GSE226353_ex_count.csv"

# File format: semicolon-separated, first column index, 2nd Geneid, 3rd Length, then samples # nolint
raw <- read.csv(count_file, sep = ";", header = TRUE, stringsAsFactors = FALSE)

# Extract gene IDs and counts (drop index and Length columns)
gene_ids <- raw$Geneid
counts <- raw[, -(1:3)]  # Samples start from column 4
rownames(counts) <- gene_ids

# Map Ensembl IDs to gene symbols
# gene_symbols <- mapIds(org.Mm.eg.db,
#                        keys = gene_ids,
#                        column = "SYMBOL",
#                        keytype = "ENSEMBL",
#                        multiVals = "first")

# # Use gene symbols as rownames (keep Ensembl ID if no symbol found)
# rownames(counts) <- ifelse(is.na(gene_symbols), gene_ids, make.unique(as.character(gene_symbols))) # nolint

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

# dds_input

############################################################
# 4. Run DESeq2
############################################################

dds <- DESeq(dds_input)

############################################################
# 5. VST transform
############################################################

vsd <- vst(dds, blind = FALSE)

############################################################
# 6. Boxplot of Reno-Risk genes (3-group)
############################################################

# Define RenoRisk genes with Ensembl IDs and symbols
renorisk_genes <- data.frame(
  ensembl_id = c("ENSMUSG00000025075", "ENSMUSG00000058952", 
                 "ENSMUSG00000026365", "ENSMUSG00000058806", 
                 "ENSMUSG00000026405"),
  symbol = c("Habp2", "Cfi", "Cfh", "Col13a1", "C4bp"),
  stringsAsFactors = FALSE
)

# Build mapping from Ensembl ID -> rowname used in dds/vsd
ensembl_to_rowname <- setNames(rownames(dds_input), gene_ids)

# Keep only RenoRisk genes present after DESeq2 prefiltering
available_rows <- ensembl_to_rowname[renorisk_genes$ensembl_id]
available_rows <- available_rows[!is.na(available_rows) & available_rows %in% rownames(dds)]

if (length(available_rows) == 0) {
  stop("None of the RenoRisk genes were found in the DESeq2 dataset after filtering.")
}

# VST expression values for requested genes
vst_mat <- assay(vsd)[available_rows, , drop = FALSE]

# Reverse lookup for labels
rowname_to_ensembl <- setNames(names(ensembl_to_rowname), ensembl_to_rowname)

plot_data_list <- list()
for (rn in rownames(vst_mat)) {
  ensembl_id <- rowname_to_ensembl[[rn]]
  gene_symbol <- renorisk_genes$symbol[match(ensembl_id, renorisk_genes$ensembl_id)]
  
  plot_data_list[[rn]] <- data.frame(
    Sample = colnames(vst_mat),
    Condition = factor(coldata$condition,
                       levels = c("WT", "CKD_veh", "CKD_RASi", "CKD_RASi_SGLT2i", "CKD_RASi_SGLT2i_MRA")),
    Expression = as.numeric(vst_mat[rn, ]),
    EnsemblID = ensembl_id,
    Gene = gene_symbol,
    stringsAsFactors = FALSE
  )
}

plot_data <- do.call(rbind, plot_data_list)
plot_data$ConditionLabel <- factor(
  as.character(plot_data$Condition),
  levels = c("WT", "CKD_veh", "CKD_RASi", "CKD_RASi_SGLT2i", "CKD_RASi_SGLT2i_MRA"),
  labels = c("WT", "Vehicle", "RASi", "RASi+SGLT2i", "RASi+SGLT2i+MRA")
)

# Wald-test pairwise contrasts for relevant comparisons
contrast_pairs <- combn(levels(coldata$condition), 2, simplify = FALSE)

pvalue_table <- data.frame()
for (rn in rownames(vst_mat)) {
  ensembl_id <- rowname_to_ensembl[[rn]]
  gene_symbol <- renorisk_genes$symbol[match(ensembl_id, renorisk_genes$ensembl_id)]
  
  for (pair in contrast_pairs) {
    res <- results(dds, contrast = c("condition", pair[1], pair[2]))
    if (rn %in% rownames(res)) {
      pvalue_table <- rbind(
        pvalue_table,
        data.frame(
          EnsemblID = ensembl_id,
          Gene = gene_symbol,
          Comparison = paste0(pair[1], " vs ", pair[2]),
          log2FC = res[rn, "log2FoldChange"],
          padj = res[rn, "padj"],
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

write.csv(pvalue_table, "RenoRisk_genes_pairwise_comparisons.csv", row.names = FALSE)
cat("Saved pairwise DESeq2 Wald-test results: RenoRisk_genes_pairwise_comparisons.csv\n")

plot_list <- list()
for (gene in unique(plot_data$Gene)) {
  gene_data <- plot_data[plot_data$Gene == gene, ]
  gene_eid <- unique(gene_data$EnsemblID)

  p <- ggplot(gene_data, aes(x = ConditionLabel, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85) +
    geom_jitter(width = 0.12, size = 2.3, alpha = 0.7) +
    scale_fill_manual(values = c(
      WT = "#00BFFF",
      CKD_veh = "#4169E1",
      CKD_RASi = "#C71585",
      CKD_RASi_SGLT2i = "#32CD32",
      CKD_RASi_SGLT2i_MRA = "#E31A1C"
    )) +
    labs(
      title = paste0(gene, " (", gene_eid, ")"),
      subtitle = "DESeq2-normalized VST expression",
      x = "",
      y = "VST expression"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      panel.grid.major.x = element_blank()
    )
  
  plot_list[[gene]] <- p
}

# Save individual plots
for (gene in names(plot_list)) {
  ggsave(paste0("RenoRisk_", gene, "_boxplot.png"), plot_list[[gene]],
         width = 8, height = 6, dpi = 300)
}
cat("Saved individual gene plots\n")

