# Install packages if needed
# install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")
# install.packages("pheatmap")
# BiocManager::install("org.Mm.eg.db")
# install.packages("ggplot2")
# install.packages("gridExtra")
# install.packages("reshape2")
# install.packages("ggsignif")

library(DESeq2)
library(apeglm)
library(pheatmap)
library(org.Mm.eg.db)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggsignif)

############################################################
# 0. Define genes data contansts
############################################################

# ENSMUSG00000026822, LCN2, ENSG00000148346, 3934_at
# ENSMUSG00000040405, HAVCR1, ENSG00000113249, 26762_at
# ENSMUSG00000046223, PLAUR, ENSG00000011422, 5329_at
# ENSMUSG00000030341 ,TNFRSF1A, ENSG00000067182, 7132_at
# ENSMUSG00000028599 ,TNFRSF1B, ENSG00000028137, 7133_at
# ENSMUSG00000032080 ,APOA4, ENSG00000110244, 337_at
# ENSMUSG00000015854 ,CD5L, ENSG00000073754, 922_at
# ENSMUSG00000020427 ,IGFBP3, ENSG00000146674, 3486_at
# EGF,        EGF,        1950,               Egf,        ENSMUSG00000028017
# FGFR4,      FGFR4,      2264,               Fgfr4,      ENSMUSG00000005320
# FGF23,      FGF23,      8074,               Fgfr4,      ENSMUSG00000000182
# Define RenoRisk genes with Ensembl IDs and symbols

renorisk_genes <- data.frame(
  ensembl_id = c("ENSMUSG00000026822",
                 "ENSMUSG00000040405",
                 "ENSMUSG00000046223",
                 "ENSMUSG00000030341",
                 "ENSMUSG00000028599",
                 "ENSMUSG00000032080",
                 "ENSMUSG00000015854",
                 "ENSMUSG00000020427",
                 "ENSMUSG00000028017",
                 "ENSMUSG00000005320",
                 "ENSMUSG00000000182",
                 "ENSMUSG00000032691"),

  symbol = c("LCN2",
              "HAVCR1",
              "PLAUR",
              "TNFRSF1A",
              "TNFRSF1B",
              "APOA4",
              "CD5L",
              "IGFBP3",
              "EGF",
              "FGFR4",
              "FGF23",
              "NLRP3"),

  stringsAsFactors = FALSE
)

############################################################
# 1. Load count matrix
############################################################

# Set working directory to the script's location
# ---- 1. Resolve script directory ----
script_dir <- if (rstudioapi::isAvailable()) {
  dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  dirname(normalizePath(sys.frames()[[1]]$ofile))
}

GEO_ACCESSION <- "GSE226353"

# Path to your file (adjust if needed)
# count_file <- paste0(GEO_ACCESSION, "_ex_count.csv")
count_file <- file.path(script_dir, paste0(GEO_ACCESSION, "_ex_count.csv"))

# count_file <- "GSE226353_ex_count.csv"
print(paste0("count_file: ", count_file))

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

dds <- DESeq(dds)

############################################################
# 5. VST transform
############################################################

vsd <- vst(dds, blind = FALSE)

############################################################
# 6. Boxplot of Reno-Risk genes (3-group)
############################################################

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

# write.csv(pvalue_table, "RenoRisk_genes_pairwise_comparisons.csv", row.names = FALSE)
# cat("Saved pairwise DESeq2 Wald-test results: RenoRisk_genes_pairwise_comparisons.csv\n")

format_pval <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("p<0.001")
  paste0("p=", formatC(round(p, 3), digits = 3, format = "f"))
}

# Comparisons: WT (reference) vs each treatment, expressed as ConditionLabel pairs
ref_comparisons <- list(
  c("WT", "Vehicle"),
  c("WT", "RASi"),
  c("WT", "RASi+SGLT2i"),
  c("WT", "RASi+SGLT2i+MRA")
)
ref_comp_orig <- c("CKD_veh", "CKD_RASi", "CKD_RASi_SGLT2i", "CKD_RASi_SGLT2i_MRA")

plot_list <- list()
for (gene in unique(plot_data$Gene)) {
  gene_data  <- plot_data[plot_data$Gene == gene, ]
  gene_eid   <- unique(gene_data$EnsemblID)
  gene_pvals <- pvalue_table[pvalue_table$Gene == gene, ]

  annotations <- vapply(ref_comp_orig, function(cond) {
    row <- gene_pvals[gene_pvals$Comparison == paste0("WT vs ", cond), ]
    if (nrow(row) > 0) format_pval(row$padj[1]) else "NS"
  }, character(1))

  y_max  <- max(gene_data$Expression, na.rm = TRUE)
  y_step <- diff(range(gene_data$Expression, na.rm = TRUE)) * 0.13
  y_pos  <- y_max + y_step * seq_along(ref_comparisons)

  p <- ggplot(gene_data, aes(x = ConditionLabel, y = Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85) +
    geom_jitter(width = 0.12, size = 2.3, alpha = 0.7) +
    geom_signif(
      comparisons = ref_comparisons,
      annotations = annotations,
      y_position  = y_pos,
      tip_length  = 0.01,
      textsize    = 3.5,
      vjust       = -0.2
    ) +
    scale_fill_manual(values = c(
      WT                  = "#00BFFF",
      CKD_veh             = "#4169E1",
      CKD_RASi            = "#C71585",
      CKD_RASi_SGLT2i     = "#32CD32",
      CKD_RASi_SGLT2i_MRA = "#E31A1C"
    )) +
    labs(
      title    = paste("Gene Expression \u2013", gene),
      subtitle = paste0("Dataset: ", GEO_ACCESSION),
      x = "Treatment",
      y = "VST-normalised expression"
    ) +
    coord_cartesian(ylim = c(NA, max(y_pos) + y_step * 0.8)) +
    theme_bw(base_size = 13) +
    theme(
      axis.text.x        = element_text(angle = 30, hjust = 1, face = "bold"),
      legend.position    = "none",
      plot.title         = element_text(face = "bold"),
      panel.grid.major.x = element_blank()
    )

  plot_list[[gene]] <- p
}

# Save individual plots
out_dir <- file.path(script_dir, "Boxplots_Other")
dir.create(out_dir, showWarnings = FALSE)
for (gene in names(plot_list)) {
  ggsave(file.path(out_dir, paste0(GEO_ACCESSION, "_boxplot_", gene, ".png")),
         plot_list[[gene]], width = 8, height = 6, dpi = 300)
}
cat("Saved individual gene plots\n")
