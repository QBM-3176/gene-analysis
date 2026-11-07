################################################################################
# GSE104954 – Kidney Tubulointerstitial Microarray Analysis
# Dataset : Affymetrix microarray, RMA-normalized, log2-scale (already processed)
# Samples : 60 human kidney biopsies across 7 disease groups
# Analysis: limma (Linear Models for Microarray Data)
#           — the correct DESeq2-analogue for continuous microarray data
# Output  : Differential expression tables + biomarker heatmap
################################################################################

# ── 0. Install / load packages ────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

for (pkg in c("limma", "GEOquery", "pheatmap", "ggplot2", "RColorBrewer")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

library(limma)
library(GEOquery)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)

# --- 0. Set working directory where the R file is located

if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
}

# --- 1. Load & parse the GEO series matrix file
GSE_ID <- "GSE104954"
GSE_ID_EXT <- "GSE104954-GPL22945"
FILE <- paste0(GSE_ID_EXT, "_series_matrix.txt")
TOP_GENES_COUNT <- 50 # 100  # Number of top variable genes for heatmap

cat("Reading series matrix …\n")
gse <- getGEO(filename = FILE, GSEMatrix = TRUE)

expr  <- exprs(gse)
rownames(expr) <- sub("_at$", "", rownames(expr))  # strip "_at" suffix
pdata <- pData(gse)

cat(sprintf("Expression matrix: %d probes × %d samples\n", nrow(expr), ncol(expr)))

# --- 2. Extract & clean sample metadata
sample_titles <- pdata$title
diagnosis <- sub("^H7-Tub-([A-Za-z&]+)[0-9]+$", "\\1", sample_titles)
diagnosis[diagnosis == "FSGS&MCD"] <- "FSGS_MCD"
diagnosis <- factor(diagnosis,
                    levels = c("LD", "DN", "FSGS", "FSGS_MCD", "MCD", "RPGN", "TN"))

cat("Sample counts per group:\n")
print(table(diagnosis))

# --- 3. Filter low-expression probes
keep <- rowSums(expr > median(expr)) >= 0.10 * ncol(expr)
expr_filt <- expr[keep, ]
cat(sprintf("After filtering: %d probes retained\n", nrow(expr_filt)))

# --- 4. Build design matrix & fit limma model
design <- model.matrix(~ 0 + diagnosis)
colnames(design) <- levels(diagnosis)
fit <- lmFit(expr_filt, design)

# --- 5. Define contrasts (Disease vs Living Donor reference)
contrast_matrix <- makeContrasts(
  DN_vs_LD       = DN       - LD,
  FSGS_vs_LD     = FSGS     - LD,
  FSGS_MCD_vs_LD = FSGS_MCD - LD,
  MCD_vs_LD      = MCD      - LD,
  RPGN_vs_LD     = RPGN     - LD,
  TN_vs_LD       = TN       - LD,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)  # empirical Bayes shrinkage — analogous to DESeq2

# --- 6. Extract differentially expressed genes (DEGs)
contrasts_list <- colnames(contrast_matrix)
de_results     <- list()

cat("\n── Differentially Expressed Genes (FDR < 0.05, |logFC| > 0.5) ──\n")
for (ctr in contrasts_list) {
  tt  <- topTable(fit2, coef = ctr, number = Inf, adjust.method = "BH", sort.by = "B")
  sig <- tt[tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.5, ]
  de_results[[ctr]] <- tt
  cat(sprintf("  %-22s : %4d DEGs  (up: %d, down: %d)\n",
              ctr, nrow(sig), sum(sig$logFC > 0), sum(sig$logFC < 0)))
}

# # --- 7. Write DE tables to CSV
# dir.create("DE_results", showWarnings = FALSE)
# for (ctr in contrasts_list) {
#   write.csv(de_results[[ctr]],
#             file = file.path("DE_results", paste0(ctr, "_limma_results.csv")),
#             row.names = TRUE)
# }
# cat("\nDE result tables written to ./DE_results/\n")

# --- 8. Select top 50 most variable genes for heatmap
gene_vars <- apply(expr_filt, 1, var)
sig_genes <- names(sort(gene_vars, decreasing = TRUE))[1:TOP_GENES_COUNT]

cat(sprintf("\nHeatmap will display %d most variable genes\n", length(sig_genes)))

# --- 9. Build heatmap expression matrix (z-score per gene)
hmap_mat <- expr_filt[sig_genes, , drop = FALSE]

# Z-score per gene, then transpose: rows = samples, cols = genes
hmap_z_genes <- t(scale(t(hmap_mat)))

diag_order  <- order(diagnosis)
diag_sorted <- diagnosis[diag_order]
hmap_z      <- t(hmap_z_genes[, diag_order])   # samples × genes

# --- 10. Annotation data frames for pheatmap
orig_sample_ids <- rownames(hmap_z)   # unique GSM IDs — capture before renaming

row_ann <- data.frame(
  Diagnosis  = diag_sorted,
  row.names  = orig_sample_ids         # unique IDs — no duplicate rowname error
)

# display labels keyed by GSM ID for pheatmap row matching
row_labels <- setNames(as.character(diag_sorted), orig_sample_ids)

diag_colours <- setNames(
  brewer.pal(nlevels(diagnosis), "Set2"),
  levels(diagnosis)
)
ann_colours <- list(Diagnosis = diag_colours)

group_sizes <- as.integer(table(diag_sorted)[levels(diagnosis)])
group_sizes  <- group_sizes[!is.na(group_sizes) & group_sizes > 0]
gaps_row     <- cumsum(head(group_sizes, -1))

# --- 11. Draw & save the heatmap
n_samples <- nrow(hmap_z)
n_genes   <- ncol(hmap_z)

png(paste0(GSE_ID_EXT, "_heatmap.png"),
    width  = max(3000, n_genes   * 55) + 500,
    height = max(3000, n_samples * 38),
    res    = 300)

pheatmap(
  hmap_z,
  scale = "column",  # Z-score per gene (column-centered)
  annotation_row    = row_ann,
  annotation_colors = ann_colours,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  labels_row        = row_labels,
  fontsize          = 15,
  fontsize_row      = 13,
  fontsize_col      = 13,
  angle_col         = 90,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  clustering_method = "complete", # "ward.D2",
  # gaps_row          = gaps_row,
  color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks            = seq(-2.5, 2.5, length.out = 101),
  border_color      = NA,
  main = paste0(
    "Biomarker Expression Profiles – ", GSE_ID_EXT, " series\n",
    "Human Kidney Tubulointerstitium (limma, z-score per gene)\n",
    n_genes, " top DEGs | FDR < 0.05 | |logFC| > 0.5 | Rows grouped by diagnosis"
  ),
  # filename = paste0(GSE_ID_EXT, "_heatmap.png"),
  # width = 10, height = 10
)

dev.off()
cat(paste0("Heatmap saved to: ", GSE_ID_EXT, "_heatmap.png\n"))

# # --- 12. Volcano plot (DN vs LD)
# volcano_df <- de_results[["DN_vs_LD"]]
# volcano_df$sig <- ifelse(
#   volcano_df$adj.P.Val < 0.05 & abs(volcano_df$logFC) > 0.5,
#   ifelse(volcano_df$logFC > 0, "Up in DN", "Down in DN"), "NS"
# )
# top_label <- head(volcano_df[order(volcano_df$adj.P.Val), ], 20)

# p_volcano <- ggplot(volcano_df, aes(x = logFC, y = -log10(adj.P.Val), colour = sig)) +
#   geom_point(alpha = 0.4, size = 1) +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
#   geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", colour = "grey40") +
#   ggrepel::geom_text_repel(
#     data = top_label, aes(label = rownames(top_label)),
#     size = 2.5, max.overlaps = 20, colour = "black"
#   ) +
#   scale_colour_manual(values = c("Up in DN" = "#d62728",
#                                  "Down in DN" = "#1f77b4", "NS" = "grey70")) +
#   labs(title = "Volcano Plot – Diabetic Nephropathy vs Living Donor",
#        subtitle = "limma moderated t-test | BH-adjusted p-value",
#        x = "log2 Fold Change", y = "-log10(FDR)", colour = NULL) +
#   theme_bw(base_size = 12) + theme(legend.position = "top")

# ggsave("GSE104954_volcano_DN_vs_LD.png", p_volcano, width = 8, height = 6, dpi = 220)
# cat("Volcano plot saved to: GSE104954_volcano_DN_vs_LD.png\n")

# # ── 13. Session info ──────────────────────────────────────────────────────────
# cat("\n── Session Info ──────────────────────────────────────────────────────\n")
# sessionInfo()