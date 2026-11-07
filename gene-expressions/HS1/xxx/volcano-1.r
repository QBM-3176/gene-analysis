# =============================================================================
# Differential Expression Analysis: Each Disease Group vs Live Donor (LD)
# Dataset : GSE104948-GPL22945 (Kidney glomerular transcriptomics)
# Method  : limma eBayes — the correct analogue of DESeq2 for log2-RMA data
#
# Groups compared (each vs LD):
#   DN      - Diabetic Nephropathy
#   FSGS    - Focal Segmental Glomerulosclerosis
#   FSGSMCD - FSGS & Minimal Change Disease (combined)
#   MCD     - Minimal Change Disease
#   RPGN    - ANCA-associated Vasculitis / Rapidly Progressive GN
#   TN      - Tumor Nephrectomy
#
# Output: one volcano PNG per comparison + one combined results CSV
# =============================================================================

# -- 0. Install / load packages -----------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

for (pkg in c("limma", "GEOquery")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
}
for (pkg in c("ggplot2", "ggrepel")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(limma)
library(GEOquery)
library(ggplot2)
library(ggrepel)

# -- 1. Load expression matrix via getGEO ------------------------------------
if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(dirname(normalizePath(sys.frames()[[1]]$ofile)))
}

GSE_ID_EXT <- "GSE104954-GPL22945"
FILE <- paste0(GSE_ID_EXT, "_series_matrix.txt")

cat("Reading GEO series matrix ...\n")
gse   <- getGEO(filename = FILE, GSEMatrix = TRUE)
expr_mat <- exprs(gse)          # probes x samples, already log2-normalized
pdata    <- pData(gse)

# Strip trailing "_at" suffix from Affymetrix probe/gene IDs
rownames(expr_mat) <- sub("_at$", "", rownames(expr_mat))

cat(sprintf("Expression matrix: %d genes x %d samples\n",
            nrow(expr_mat), ncol(expr_mat)))

# -- 2. Derive group labels from sample titles ---------------------------------
# Titles look like: H7-Glom-RPGN902, H7-Glom-LD_974, H7-Glom-FSGS&MCD905
# Strategy: strip trailing digits, take token after last "-", remove non-alpha
sample_titles <- pdata$title

get_group <- function(title) {
  core <- sub("[0-9]+$", "", title)      # drop trailing digits
  grp  <- sub(".*-", "", core)           # keep text after last "-"
  grp  <- gsub("[^A-Za-z]", "", grp)    # remove "&", "_" etc. -> FSGSMCD
  return(grp)
}

groups <- sapply(sample_titles, get_group)
names(groups) <- colnames(expr_mat)

cat("\nSample counts per group:\n")
print(table(groups))

disease_groups <- sort(setdiff(unique(groups), "LD"))
cat(sprintf("\nComparisons to run: %s vs LD\n\n",
            paste(disease_groups, collapse = ", ")))

# -- 3. Low-variance filter (applied once, globally) --------------------------
row_vars  <- apply(expr_mat, 1, var, na.rm = TRUE)
expr_filt <- expr_mat[row_vars > quantile(row_vars, 0.10, na.rm = TRUE), ]
cat(sprintf("After low-variance filter: %d genes retained\n\n", nrow(expr_filt)))

# -- 4. Volcano plot helper ---------------------------------------------------
make_volcano <- function(res_df, grp_label, n_disease, n_ld,
                         fc_thresh = 1.0, fdr_thresh = 0.05) {

  res_df$neg_log10_p <- -log10(res_df$pvalue)

  sig_up   <- sum(res_df$padj < fdr_thresh & res_df$log2FC >  fc_thresh, na.rm = TRUE)
  sig_down <- sum(res_df$padj < fdr_thresh & res_df$log2FC < -fc_thresh, na.rm = TRUE)

  res_df$significance <- "Not significant"
  res_df$significance[res_df$padj < fdr_thresh & res_df$log2FC >  fc_thresh] <- "Up"
  res_df$significance[res_df$padj < fdr_thresh & res_df$log2FC < -fc_thresh] <- "Down"
  res_df$significance <- factor(res_df$significance,
                                levels = c("Up", "Down", "Not significant"))

  # Top 15 DEGs to label (joint rank on p-value and |FC|)
  res_df$rank_score <- rank(res_df$pvalue) + rank(-abs(res_df$log2FC))
  label_df <- head(res_df[order(res_df$rank_score), ], 15)

  sig_colours <- c(
    "Up"              = "#E63946",
    "Down"            = "#457B9D",
    "Not significant" = "#ADB5BD"
  )

  x_lim <- max(abs(res_df$log2FC), na.rm = TRUE) * 1.08
  y_max <- max(res_df$neg_log10_p, na.rm = TRUE)

  ggplot(res_df, aes(x = log2FC, y = neg_log10_p,
                     colour = significance, size = significance)) +
    geom_point(alpha = 0.65, shape = 16) +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh),
               linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = -log10(fdr_thresh),
               linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_text_repel(
      data           = label_df,
      aes(label      = GeneID),
      size           = 2.7,
      fontface       = "bold",
      colour         = "black",
      box.padding    = 0.4,
      max.overlaps   = 20,
      segment.colour = "grey50",
      segment.size   = 0.3
    ) +
    scale_colour_manual(values = sig_colours) +
    scale_size_manual(values  = c("Up" = 1.8, "Down" = 1.8,
                                  "Not significant" = 0.8)) +
    xlim(-x_lim, x_lim) +
    annotate("text", x =  x_lim * 0.72, y = y_max * 0.97,
             label = paste0("Up: ", sig_up),
             colour = "#E63946", fontface = "bold", size = 4) +
    annotate("text", x = -x_lim * 0.72, y = y_max * 0.97,
             label = paste0("Down: ", sig_down),
             colour = "#457B9D", fontface = "bold", size = 4) +
    labs(
      title    = sprintf("DEG: %s vs Live Donor (LD)", grp_label),
      subtitle = sprintf(
        "limma eBayes | BH FDR < %g & |log2FC| > %g  |  n(%s) = %d,  n(LD) = %d",
        fdr_thresh, fc_thresh, grp_label, n_disease, n_ld),
      x        = paste0("log2 Fold Change (", grp_label, " / LD)"),
      y        = "-log10(p-value)",
      colour   = "Significance",
      caption  = "Labels: Entrez Gene IDs of top 15 most significant DEGs"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(size = 9.5, colour = "grey30"),
      legend.position  = "right",
      legend.title     = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    guides(size = "none")
}

# -- 5. Loop: one limma analysis per disease group ----------------------------
all_results <- list()

for (grp in disease_groups) {

  cat(sprintf("== %s vs LD ==\n", grp))

  # Subset to this group + LD only
  keep      <- groups %in% c(grp, "LD")
  sub_expr  <- expr_filt[, keep]
  sub_group <- groups[keep]

  n_disease <- sum(sub_group == grp)
  n_ld      <- sum(sub_group == "LD")
  cat(sprintf("   n(%s) = %d   n(LD) = %d\n", grp, n_disease, n_ld))

  # Design matrix with LD as the reference level
  sub_factor <- factor(sub_group, levels = c("LD", grp))
  design     <- model.matrix(~ sub_factor)
  colnames(design) <- c("Intercept", paste0(grp, "_vs_LD"))

  # Fit limma (analogous to DESeq2: model fit + dispersion + Wald test)
  fit  <- lmFit(sub_expr, design)
  fit2 <- eBayes(fit)   # empirical Bayes moderated t-statistics

  # Full results table with BH FDR correction
  res <- topTable(
    fit2,
    coef          = paste0(grp, "_vs_LD"),
    number        = Inf,
    adjust.method = "BH",
    sort.by       = "P"
  )

  # Rename columns to match DESeq2 conventions
  colnames(res)[colnames(res) == "logFC"]     <- "log2FC"
  colnames(res)[colnames(res) == "P.Value"]   <- "pvalue"
  colnames(res)[colnames(res) == "adj.P.Val"] <- "padj"
  res$GeneID     <- rownames(res)
  res$comparison <- paste0(grp, "_vs_LD")

  sig_up   <- sum(res$padj < 0.05 & res$log2FC >  1, na.rm = TRUE)
  sig_down <- sum(res$padj < 0.05 & res$log2FC < -1, na.rm = TRUE)
  cat(sprintf("   Significant (padj<0.05 & |log2FC|>1): %d up, %d down\n\n",
              sig_up, sig_down))

  all_results[[grp]] <- res

  # Build and save volcano
  p <- make_volcano(res, grp, n_disease, n_ld)
  out_file <- paste0(GSE_ID_EXT, "_volcano_", grp, "_vs_LD.png")
  ggsave(out_file, plot = p, width = 10, height = 7, dpi = 300, bg = "white")
  cat(sprintf("   Volcano saved -> %s\n\n", out_file))
}

# -- 6. Save combined results CSV ---------------------------------------------
# combined <- do.call(rbind, lapply(all_results, function(df) {
#   df[, c("comparison", "GeneID", "log2FC", "AveExpr", "t", "pvalue", "padj", "B")]
# }))
# rownames(combined) <- NULL

# write.csv(combined, "limma_results_all_groups_vs_LD.csv", row.names = FALSE)
# cat("Combined results saved -> limma_results_all_groups_vs_LD.csv\n")
# cat("All done.\n")
