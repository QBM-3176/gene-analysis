############################################################
# limma pipeline + Heatmap for GSE104954-GPL22945
# Input: GSE104954-GPL22945_series_matrix.txt
#        (Affymetrix microarray, RMA-normalized log2 — pre-processed)
# Analysis: limma (correct DESeq2-analogue for microarray data)
# Output: Results CSVs + heatmap "GSE104954-GPL22945_heatmap_topXX.png"
############################################################

# Install packages if needed
# install.packages("BiocManager")
# BiocManager::install("limma")
# BiocManager::install("GEOquery")
# install.packages("pheatmap")

library(limma)
library(GEOquery)
library(pheatmap)

############################################################
# 1. Load expression matrix
############################################################

if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
}

GSE_ID_EXT <- "GSE104954-GPL22945"
FILE <- paste0(GSE_ID_EXT, "_series_matrix.txt")

TOP_GENES_COUNT <- 50  # Number of top variable genes for heatmap

cat("Reading series matrix …\n")
gse   <- getGEO(filename = FILE, GSEMatrix = TRUE)
expr  <- exprs(gse)          # probes × samples, already log2-normalized
pdata <- pData(gse)

cat(sprintf("Expression matrix: %d probes × %d samples\n", nrow(expr), ncol(expr)))

############################################################
# 2. Extract sample metadata (colData equivalent)
############################################################

sample_titles <- pdata$title
diagnosis <- sub("^H7-Tub-([A-Za-z&]+)[0-9]+$", "\\1", sample_titles)
diagnosis[diagnosis == "FSGS&MCD"] <- "FSGS_MCD"
diagnosis <- factor(diagnosis,
                    levels = c("LD", "DN", "FSGS", "FSGS_MCD", "MCD", "RPGN", "TN"))

cat("Sample counts per group:\n")
print(table(diagnosis))

############################################################
# 3. Filter low-expression probes (analogous to DESeq2 prefilter)
############################################################

keep   <- rowSums(expr > median(expr)) >= 0.10 * ncol(expr)
expr_f <- expr[keep, ]
cat(sprintf("After filtering: %d probes retained\n", nrow(expr_f)))

############################################################
# 4. Build design matrix & fit limma model (analogous to DESeqDataSet + DESeq())
############################################################

design <- model.matrix(~ 0 + diagnosis)
colnames(design) <- levels(diagnosis)
fit <- lmFit(expr_f, design)

############################################################
# 5. Define contrasts (Disease vs Living Donor — analogous to results())
############################################################

# contrast_matrix <- makeContrasts(
#   DN_vs_LD       = DN       - LD,
#   FSGS_vs_LD     = FSGS     - LD,
#   FSGS_MCD_vs_LD = FSGS_MCD - LD,
#   MCD_vs_LD      = MCD      - LD,
#   RPGN_vs_LD     = RPGN     - LD,
#   TN_vs_LD       = TN       - LD,
#   levels = design
# )

# fit2 <- contrasts.fit(fit, contrast_matrix)
# fit2 <- eBayes(fit2)   # empirical Bayes shrinkage — analogous to DESeq2 shrinkage

############################################################
# 6. Extract differentially expressed genes (analogous to results() + write.csv)
############################################################

# contrasts_list <- colnames(contrast_matrix)

# cat("\n── Differentially Expressed Genes (FDR < 0.05, |logFC| > 0.5) ──\n")
# for (ctr in contrasts_list) {
#   tt  <- topTable(fit2, coef = ctr, number = Inf, adjust.method = "BH", sort.by = "B")
#   sig <- tt[tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.5, ]
#   write.csv(tt, file = paste0("limma_", ctr, ".csv"))
#   cat(sprintf("  %-22s : %4d DEGs  (up: %d, down: %d)\n",
#               ctr, nrow(sig), sum(sig$logFC > 0), sum(sig$logFC < 0)))

############################################################
# 7. (No VST needed — microarray data is already log2-normalized)
############################################################

# expr_f is already in log2 space — used directly for the heatmap

############################################################
# 8. Heatmap: Top TOP_GENES_COUNT most variable genes
############################################################

all_genes_sorted <- names(sort(apply(expr_f, 1, var), decreasing = TRUE))
top_var_genes    <- all_genes_sorted[1:TOP_GENES_COUNT]

# Subset expression matrix
mat <- expr_f[top_var_genes, ]

# Transpose matrix (genes become columns, samples become rows)
mat <- t(mat)

############################################################
# 9. Annotation (analogous to annotation_row in RDaseq2.r)
############################################################

annotation_row <- data.frame(
  Diagnosis = diagnosis,
  row.names = rownames(mat)
)

# Map sample IDs → diagnosis labels for row display
row_labels <- setNames(as.character(diagnosis), rownames(mat))

annotation_colors <- list(
  Diagnosis = c(
    LD       = "#00FFFF",
    DN       = "#0000FF",
    FSGS     = "#FF00FF",
    FSGS_MCD = "#00FF00",
    MCD      = "#FFA500",
    RPGN     = "#FF0000",
    TN       = "#FFFF00"
  )
)

############################################################
# 10. Draw & save heatmap (same style as RDaseq2.r)
############################################################

out_file  <- paste0(GSE_ID_EXT, "_heatmap_top", TOP_GENES_COUNT, ".png")
cell_size <- 10  # points — same for width & height → square cells

# Total width: cells + row dendrogram + row labels + annotation + legend
# pheatmap honours cellwidth/cellheight and auto-sizes internally;
# width/height are the canvas upper bound (in inches at 300 dpi).
# 22 × 300 = 6600 px — ample room for legend text "FSGS_MCD".
pheatmap(mat,
         scale                    = "column",  # Z-score per gene
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method        = "complete",
         cellwidth  = cell_size,
         cellheight = cell_size,
         color = colorRampPalette(c("navy", "white", "red"))(TOP_GENES_COUNT),
         breaks = seq(-1.5, 1.5, length.out = TOP_GENES_COUNT + 1),
         annotation_row    = annotation_row,
         annotation_colors = annotation_colors,
         show_rownames = TRUE,
         labels_row    = row_labels,
         show_colnames = TRUE,
         fontsize      = 14,   # title = fontsize + 4 = 20 pt
         fontsize_row  = 11,
         fontsize_col  = 11,
         main = paste0(GSE_ID_EXT, "\nTop ", TOP_GENES_COUNT,
                       " Variable Genes (Limma analysis)\n"),
         filename = out_file,
         width = 12, height = 12)

cat(paste0("Heatmap saved as ", out_file, "\n"))
