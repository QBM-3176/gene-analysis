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

library(DESeq2)
library(apeglm)
library(pheatmap)
library(org.Mm.eg.db)

############################################################
# 1. Load count matrix
############################################################

# Path to your file (adjust if needed)
count_file <- "GSE226353_ex_count.csv"

TOP_GENES_COUNT <- 50 # 100  # Number of top variable genes for heatmap

# File format: semicolon-separated, first column index, 2nd Geneid, 3rd Length, then samples # nolint
raw <- read.csv(count_file, sep = ";", header = TRUE, stringsAsFactors = FALSE)

# Extract gene IDs and counts (drop index and Length columns)
gene_ids <- raw$Geneid
counts <- raw[, -(1:3)]  # Samples start from column 4
rownames(counts) <- gene_ids

# Map Ensembl IDs to gene symbols
gene_symbols <- mapIds(org.Mm.eg.db,
                       keys = gene_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Use gene symbols as rownames (keep Ensembl ID if no symbol found)
rownames(counts) <- ifelse(is.na(gene_symbols), gene_ids, make.unique(as.character(gene_symbols))) # nolint

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
# 5. Extract shrunken results (all treatments vs WT)
############################################################

contrasts <- c("CKD_veh", "CKD_RASi", "CKD_RASi_SGLT2i", "CKD_RASi_SGLT2i_MRA")

for (ct in contrasts) {
  res <- results(dds, contrast = c("condition", ct, "WT"))
  res_shrunk <- lfcShrink(dds, coef = paste0("condition_", ct, "_vs_WT"), type = "apeglm")
  res_shrunk <- res_shrunk[order(res_shrunk$padj), ]
  write.csv(as.data.frame(res_shrunk), file = paste0("DESeq2_", ct, "_vs_WT.csv"))
  cat("Saved", paste0("DESeq2_", ct, "_vs_WT.csv"), "\n")
}

############################################################
# 6. QC plots
############################################################

plotDispEsts(dds)
plotMA(results(dds, contrast = c("condition", "CKD_RASi_SGLT2i_MRA", "WT")), ylim = c(-8, 8))

############################################################
# 7. VST transform for heatmap
############################################################

vsd <- vst(dds, blind = FALSE)

############################################################
# 8. Heatmap: Top TOP_GENES_COUNT most variable genes
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
pheatmap(mat,
         scale = "column",  # Z-score per gene (column-centered)
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("navy", "white", "red"))(TOP_GENES_COUNT),
         breaks = seq(-1.5, 1.5, length.out = TOP_GENES_COUNT+1),  # Exact -1.5 to +1.5 range
         annotation_row = annotation_row,
         annotation_colors = annotation_colors,
         show_rownames = TRUE,  # Show sample labels
         show_colnames = TRUE,  # Show gene labels
         fontsize_col = 8,
         main = paste0("GSE226353: Top ",  TOP_GENES_COUNT, " Variable Genes"),
         filename = paste0("GSE226353_heatmap_top", TOP_GENES_COUNT, ".png"),  # Saves PNG
         width = 12, height = 4)

cat(paste0("Heatmap saved as GSE226353_heatmap_top", TOP_GENES_COUNT, ".png\n"))

############################################################
# 9. Normalized counts
############################################################

norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, file = "GSE226353_normalized_counts.csv")
