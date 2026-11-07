#!/usr/bin/env Rscript
# =============================================================================
# MetaAnalysis — pooling strategy (plan §5): random-effects meta-analysis of the 5-gene panel
#
# SOURCE: CKD_Biomarker_Pooled_Analysis_Plan_Final.pdf — "PROJECT ANALYSIS PLAN &
#   TECHNICAL MASTER PLAN", RenoRisk Validation Across Multiple CKD Indications,
#   Version 2.0 (Final Analysis Plan), July 2026. Project Lead: Amela Jusic;
#   Execution Lead: Stojanko Tudjarski. Quotations below are verbatim.
#
# PLAN REQUIREMENTS IMPLEMENTED HERE
#   §5 Methodology: Pooling Strategy for C4BPA, COL13A1, HABP2, CFH, CFI:
#     "To synthesize findings across Microarray and RNA-Seq cohorts without
#     platform bias, we apply summary-level random-effects meta-analysis as the
#     primary framework, complemented by mega-analysis within platform tiers."
#     Random-Effects Meta-Analysis Model:
#     "Per-Study Effect Size: For each RenoRisk gene (C4BPA, COL13A1, HABP2, CFH,
#     CFI), compute Hedges' g standardized mean difference (or log2 fold change)
#     and variance in each dataset."
#     "Pooled Summary Estimate: Combine study effect sizes using DerSimonian-Laird
#     random-effects models via the metafor package in R."
#     "Heterogeneity Evaluation: Evaluate Cochran's Q and I² statistics to measure
#     biological and technical heterogeneity across cohorts."
#   §7 Plot 2 — Forest Plots for 5 RenoRisk Genes:
#     "Individual study effect sizes and pooled summary diamonds for C4BPA,
#     COL13A1, HABP2, CFH, and CFI."
#     Example, Annex 2 Figure 1 ("Forest plot example"): "Forest Plot Example:
#     COL13A1 Across Key Cohorts", one row per cohort plus "Pooled Diamond (RE
#     Model)", x axis "Log2 Fold Change (COL13A1: Disease vs. Control)".
#   §7 Plot 4 — Indication Stratified Expression Heatmaps:
#     "Expression of C4BPA, COL13A1, HABP2, CFH, CFI across CKD, DKD, IgAN, FSGS,
#     AKI, and Healthy Controls."
#     Example, Annex 2 Figure 3 ("Heat map – expression profile of RenoRisk across
#     priority indications"): genes x {Control, CKD Stg 1-2, CKD Stg 3-5, DKD,
#     IgAN, FSGS, AKI}, coloured by "Z-Score Expression".
#   §8 R code template ("Random-Effects Meta-Analysis for COL13A1 across cohorts"):
#     meta_col13a1 <- rma(yi = col13a1_yi, vi = col13a1_vi, method = "DL")
#     summary(meta_col13a1)
#     forest(meta_col13a1, slab = c(...), main = "COL13A1 Expression Across
#     Kidney Disease Cohorts")
#   §2 Target RenoRisk Biomarker Panel Specification:
#     "The RenoRisk panel consists of 5 specific target biomarkers categorized into
#     primary diagnostic candidates and additional complement system regulators"
#     Main RenoRisk: C4BPA (C4b-Binding Protein Alpha), COL13A1 (Collagen Type XIII
#     Alpha 1 Chain), HABP2 (Hyaluronan Binding Protein 2); Additional RenoRisk:
#     CFH (Complement Factor H), CFI (Complement Factor I).
#
# HOW THIS SCRIPT ADDRESSES THEM (and where it deviates, with the reason)
#   * Per-study effect: Hedges' g (escalc "SMD", bias-corrected) is primary and
#     log2 fold change (escalc "MD") secondary — both named by the plan. g is
#     primary because it is scale-free: four cohorts are per-gene centred, so raw
#     log2 differences are not on a common footing.
#   * "each dataset": the unit is the Stage 2/3 batch unit (15 microarray + 3 RNA-Seq = 18, all pooled; GSE30122 split
#     into its glomerular and tubular parts) because a dataset holding two
#     compartments has two different disease effects.
#   * Pooling: metafor::rma(method = "DL") as in the plan and its §8 template —
#     overall, per compartment, and a compartment moderator test.
#   * Heterogeneity: Cochran's Q (and its p), I², plus tau² and 95% prediction
#     intervals; leave-one-out sensitivity on the overall model.
#   * Microarray vs RNA-Seq "without platform bias": RNA-Seq studies enter as
#     per-study effect sizes only (summary level; g is scale-free), never as
#     pooled sample values. One RNA-Seq series has a Disease-vs-Control contrast:
#     GSE142025 (whole biopsy; DKD vs tumor-nephrectomy controls; C4BPA not in
#     its deposit). It is a regular Stage 2/3 cohort (panel_uncorrected_long).
#   * Plot 2: one forest plot per gene (study rows, subgroup diamonds, overall RE
#     diamond), on Hedges' g rather than the example's log2FC (reason above);
#     plus a one-page summary forest of all genes x compartments.
#   * Plot 4: rows are the indications present in the data (DKD, FSGS, IgAN, MN,
#     MCD, ANCA vasculitis, lupus nephritis, ... and CKD unspecified); AKI is not
#     in any curated cohort and no cohort records CKD stage, so "AKI" and "CKD Stg
#     1-2 / 3-5" cannot be shown. Cells are within-study Hedges' g of each
#     indication vs that study's own controls, pooled by DL — not z-scored raw
#     expression, because absolute levels are not comparable across cohorts
#     (per-gene-centred deposits, different platforms).
#   * Not in the plan, added as a sensitivity check: whether the pooled effect
#     depends on the control tissue (living donor vs tumor nephrectomy).
#
# Plan: per-study Hedges' g (or log2 fold change) and its variance for each
# RenoRisk gene; pool with DerSimonian–Laird random-effects models (metafor);
# Cochran's Q and I² for heterogeneity. Forest plots per gene (plan Plot 2).
#
# Input : Results/Stage23/panel_uncorrected_long.xlsx (written by Stage 2/3).
#         Effect sizes are computed WITHIN each study, so batch correction is
#         neither needed nor wanted — the uncorrected values are the right input.
# Output: Results/PoolingStrategy/ — tables (.xlsx) and figures (forest plots, heatmap, sensitivity)
#
# DESIGN
#   * Study unit = Stage 2/3 batch unit: each study, with GSE30122 split into
#     its glomerular and tubular parts (15 microarray + 3 RNA-Seq units, all pooled). A unit enters a gene's
#     meta-analysis when both arms have >= cfg$min_n samples.
#   * Primary measure: Hedges' g (bias-corrected SMD, Disease minus Control).
#     It is scale-free, which matters here: cohorts were deposited on different
#     scales (four are per-gene centred), so raw log2 differences are not on a
#     common footing. log2 fold change (raw mean difference) is reported as the
#     secondary measure.
#   * Pooled three ways per gene: all units (the plan's overall estimate), per
#     compartment (subgroups, each with its own tau²), and a moderator test of
#     whether the effect differs between compartments (QM).
#   * Leave-one-out sensitivity on the overall model.
#
# CAVEATS
#   * A dataset fully contained in another (cfg$subsets; none since GSE30528/
#     GSE30529 were dropped on 2026-11-03) would get a per-study effect
#     (`Pooled` = FALSE) but stay out of every pooled estimate.
#   * Few units per subgroup (Tubulointerstitial 4, Whole biopsy 4): tau² and
#     I² there are very imprecise; DL can return tau² = 0 by default.
#   * Disease is the binary Control/Disease collapse of Stage 2/3.
#   * Values are submitter-processed (five pipelines).
# =============================================================================

cfg <- list(
  project_dir = path.expand("~/Striki/Quant/quant.2/MetaAnalysis"),
  genes       = c("C4BPA", "COL13A1", "HABP2", "CFH", "CFI"),
# GSE30528/GSE30529 (subsets of GSE30122) were dropped and renamed -dup on 2026-11-03 (by instruction); the guard stays for any future contained dataset.
  subsets     = character(0),   # datasets fully contained in another: not pooled
  min_n       = 3,       # per arm, for a unit to enter a gene's meta-analysis
  method      = "DL",    # DerSimonian–Laird, per the plan
  width       = 9,
  dpi         = 300
)
cfg$results_dir <- file.path(cfg$project_dir, "Results")   # everything the pipeline generates
cfg$in_file <- file.path(cfg$results_dir, "Stage23", "panel_uncorrected_long.xlsx")
cfg$out_dir <- file.path(cfg$results_dir, "PoolingStrategy")
cfg$fig_dir <- cfg$out_dir   # figures sit with the tables

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(metafor)
  library(readxl); library(writexl)
})
set.seed(1)
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg$fig_dir, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

if (!file.exists(cfg$in_file))
  stop("missing ", cfg$in_file, " - run Scripts/stage23_panel_boxplots.R first")

# ---------------------------------------------------------------------------
# 1. Per-study summary statistics -> effect sizes
# ---------------------------------------------------------------------------
# Column types are fixed, not guessed: GSM, StudyID, Condition, Compartment,
# Batch, Gene (text), Expression (numeric).
long <- as.data.table(read_excel(cfg$in_file, col_types = c(rep("text", 6), "numeric")))
stats <- long[, .(m = mean(Expression), s = sd(Expression), n = .N),
              by = .(Gene, Batch, StudyID, Compartment, Condition)]
wide  <- dcast(stats, Gene + Batch + StudyID + Compartment ~ Condition,
               value.var = c("m", "s", "n"))
wide  <- wide[n_Control >= cfg$min_n & n_Disease >= cfg$min_n &
              s_Control > 0 & s_Disease > 0]

# Hedges' g: escalc's SMD applies the small-sample bias correction.
es_g  <- escalc(measure = "SMD", m1i = m_Disease, sd1i = s_Disease, n1i = n_Disease,
                m2i = m_Control, sd2i = s_Control, n2i = n_Control, data = wide)
es_fc <- escalc(measure = "MD",  m1i = m_Disease, sd1i = s_Disease, n1i = n_Disease,
                m2i = m_Control, sd2i = s_Control, n2i = n_Control, data = wide)
eff <- data.table(wide[, .(Gene, Batch, StudyID, Compartment,
                           n_Control, n_Disease)],
                  g = as.numeric(es_g$yi), g_var = as.numeric(es_g$vi),
                  log2FC = as.numeric(es_fc$yi), log2FC_var = as.numeric(es_fc$vi))
eff[, Gene := factor(Gene, levels = intersect(cfg$genes, unique(Gene)))]
eff[, Pooled := !StudyID %in% cfg$subsets]   # subsets of GSE30122: not pooled
setorder(eff, Gene, Compartment, Batch)
write_xlsx(eff, file.path(cfg$out_dir, "pooling_study_effects.xlsx"))
msg(nrow(eff), " gene x unit effect sizes from ", uniqueN(eff$Batch), " units; pooled: ",
    uniqueN(eff[(Pooled), Batch]), " units")

# ---------------------------------------------------------------------------
# 2. Random-effects pooling: overall, per compartment, moderator test
# ---------------------------------------------------------------------------
summarise_rma <- function(fit, gene, measure, subset) {
  pr <- predict(fit)
  data.table(Gene = gene, Measure = measure, Subset = subset, k = fit$k,
             estimate = as.numeric(fit$beta), SE = fit$se,
             ci_lb = fit$ci.lb, ci_ub = fit$ci.ub,
             pi_lb = pr$pi.lb, pi_ub = pr$pi.ub,
             z = fit$zval, p = fit$pval,
             Q = fit$QE, Q_df = fit$k - 1, Q_p = fit$QEp,
             I2 = fit$I2, tau2 = fit$tau2)
}

pooled <- list(); moder <- list(); fits <- list()
for (gene in levels(eff$Gene)) for (ms in c("g", "log2FC")) {
  d  <- eff[Gene == gene & Pooled]
  yi <- d[[ms]]; vi <- d[[paste0(ms, "_var")]]
  fit <- rma(yi, vi, method = cfg$method)
  if (ms == "g") fits[[gene]] <- list(fit = fit, d = d)
  pooled[[length(pooled) + 1]] <- summarise_rma(fit, gene, ms, "All units")
  for (cp in sort(unique(d$Compartment))) {
    idx <- d$Compartment == cp
    pooled[[length(pooled) + 1]] <- if (sum(idx) >= 2)
      summarise_rma(rma(yi[idx], vi[idx], method = cfg$method), gene, ms, cp)
    else  # one unit: no pooling, report the study estimate itself
      data.table(Gene = gene, Measure = ms, Subset = cp, k = 1L,
                 estimate = yi[idx], SE = sqrt(vi[idx]),
                 ci_lb = yi[idx] - qnorm(0.975) * sqrt(vi[idx]),
                 ci_ub = yi[idx] + qnorm(0.975) * sqrt(vi[idx]),
                 z = yi[idx] / sqrt(vi[idx]), p = 2 * pnorm(-abs(yi[idx] / sqrt(vi[idx]))))
  }
  if (uniqueN(d$Compartment) > 1) {
    mf <- rma(yi, vi, mods = ~ Compartment, data = d, method = cfg$method)
    moder[[length(moder) + 1]] <- data.table(Gene = gene, Measure = ms, k = mf$k,
      QM = mf$QM, QM_df = mf$QMdf[1], QM_p = mf$QMp, residual_I2 = mf$I2)
  }
}
pooled <- rbindlist(pooled, fill = TRUE)
# BH across the 5 genes, within each measure x subset.
pooled[, FDR := p.adjust(p, "BH"), by = .(Measure, Subset)]
moder  <- rbindlist(moder)
moder[, FDR := p.adjust(QM_p, "BH"), by = Measure]

num <- names(pooled)[vapply(pooled, is.numeric, TRUE) & !names(pooled) %in% c("k", "Q_df")]
out_pooled <- copy(pooled)[, (num) := lapply(.SD, signif, 4), .SDcols = num]
write_xlsx(out_pooled, file.path(cfg$out_dir, "pooling_estimates.xlsx"))
write_xlsx(moder[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)],
           file.path(cfg$out_dir, "pooling_compartment_moderator.xlsx"))

msg("pooled Hedges' g, all units (DL random effects):")
print(out_pooled[Measure == "g" & Subset == "All units",
                 .(Gene, k, g = estimate, ci_lb, ci_ub, p, FDR, I2 = round(I2, 1),
                   tau2, Q_p)])
msg("pooled Hedges' g by compartment:")
print(out_pooled[Measure == "g" & Subset != "All units",
                 .(Gene, Subset, k, g = estimate, ci_lb, ci_ub, FDR, I2 = round(I2, 1))])
msg("does the effect differ between compartments? (moderator test, Hedges' g)")
print(moder[Measure == "g", .(Gene, k, QM = round(QM, 2), QM_df, QM_p = signif(QM_p, 3),
                              FDR = signif(FDR, 3), residual_I2 = round(residual_I2, 1))])

# ---------------------------------------------------------------------------
# 3. Leave-one-out sensitivity (overall model, Hedges' g)
# ---------------------------------------------------------------------------
loo <- rbindlist(lapply(names(fits), function(gene) {
  l <- leave1out(fits[[gene]]$fit)
  data.table(Gene = gene, Omitted = fits[[gene]]$d$Batch,
             g = l$estimate, ci_lb = l$ci.lb, ci_ub = l$ci.ub, p = l$pval,
             I2 = l$I2, tau2 = l$tau2)
}))
write_xlsx(loo[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)],
           file.path(cfg$out_dir, "pooling_leave1out.xlsx"))
loo_rng <- loo[, .(g_min = round(min(g), 3), g_max = round(max(g), 3),
                   max_p = signif(max(p), 3)), by = Gene]
msg("leave-one-out range of the pooled g (max p = least significant omission):")
print(loo_rng)

# ---------------------------------------------------------------------------
# 4. Forest plots — one per gene: units grouped by compartment, subgroup
#    diamonds, overall diamond
# ---------------------------------------------------------------------------
diamond <- function(est, lb, ub, y, h = 0.35)
  data.table(x = c(lb, est, ub, est), y = c(y, y + h, y, y - h))

for (gene in levels(eff$Gene)) {
  d  <- eff[Gene == gene & Pooled]           # the plot shows what is pooled
  d[, `:=`(lo = g - qnorm(0.975) * sqrt(g_var), hi = g + qnorm(0.975) * sqrt(g_var))]
  rows <- list(); dias <- list(); y <- 0
  pg <- pooled[Gene == gene & Measure == "g"]
  for (cp in rev(sort(unique(d$Compartment)))) {
    sg <- pg[Subset == cp]
    y <- y + 1
    dias[[length(dias) + 1]] <- cbind(diamond(sg$estimate, sg$ci_lb, sg$ci_ub, y),
                                      grp = cp, kind = "Subgroup")
    rows[[length(rows) + 1]] <- data.table(y = y, label = sprintf("%s (k=%d)", cp, sg$k),
      est = sg$estimate, lo = sg$ci_lb, hi = sg$ci_ub, kind = "Subgroup",
      n = "", txt = sprintf("%.2f [%.2f, %.2f]", sg$estimate, sg$ci_lb, sg$ci_ub))
    dd <- d[Compartment == cp][order(-g)]
    for (i in seq_len(nrow(dd))) {
      y <- y + 1
      rows[[length(rows) + 1]] <- data.table(y = y, label = paste0("   ", dd$Batch[i]),
        est = dd$g[i], lo = dd$lo[i], hi = dd$hi[i], kind = "Study",
        n = sprintf("%d / %d", dd$n_Control[i], dd$n_Disease[i]),
        txt = sprintf("%.2f [%.2f, %.2f]", dd$g[i], dd$lo[i], dd$hi[i]))
    }
    y <- y + 0.6
  }
  ov <- pg[Subset == "All units"]
  y <- y + 1
  dias[[length(dias) + 1]] <- cbind(diamond(ov$estimate, ov$ci_lb, ov$ci_ub, y, 0.45),
                                    grp = "All", kind = "Overall")
  rows[[length(rows) + 1]] <- data.table(y = y, label = sprintf("Overall RE model (k=%d)", ov$k),
    est = ov$estimate, lo = ov$ci_lb, hi = ov$ci_ub, kind = "Overall", n = "",
    txt = sprintf("%.2f [%.2f, %.2f]", ov$estimate, ov$ci_lb, ov$ci_ub))
  rows <- rbindlist(rows); dias <- rbindlist(dias)
  ymax <- max(rows$y)
  rows[, y := ymax + 1 - y]; dias[, y := ymax + 1 - y]   # overall on top

  xr <- range(c(rows$lo, rows$hi, 0), na.rm = TRUE)
  st <- rows[kind == "Study"]
  p <- ggplot() +
    geom_vline(xintercept = 0, colour = "grey55", linetype = 2) +
    geom_vline(xintercept = ov$estimate, colour = "#1B7F79", linewidth = 0.3, alpha = 0.6) +
    geom_errorbar(data = st, aes(y = y, xmin = lo, xmax = hi), width = 0.25,
                  orientation = "y", colour = "#1F3B5C") +
    geom_point(data = st, aes(est, y), shape = 15, size = 2.2, colour = "#1F3B5C") +
    geom_polygon(data = dias, aes(x, y, group = interaction(grp, kind), fill = kind)) +
    scale_fill_manual(values = c(Subgroup = "#7FB8B4", Overall = "#1B7F79"), guide = "none") +
    scale_y_continuous(breaks = rows$y, labels = rows$label,
                       sec.axis = dup_axis(labels = paste(rows$n, rows$txt, sep = "   "),
                                           name = "n Ctrl / Dis     g [95% CI]")) +
    labs(title = sprintf("%s — Disease vs Control, random-effects meta-analysis", gene),
         subtitle = sprintf("Hedges' g, DerSimonian–Laird.  Overall I² = %.0f%%, tau² = %.3f, Q p = %s",
                            ov$I2, ov$tau2, format.pval(ov$Q_p, digits = 2)),
         x = "Hedges' g (Disease − Control)", y = NULL,
         caption = "Units: studies (GSE30122 split by compartment). RNA-Seq: GSE142025, GSE141295, GSE175759; the rest microarrays.") +
    coord_cartesian(xlim = xr) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          axis.ticks.y = element_blank(), plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey35", size = 9),
          plot.caption = element_text(colour = "grey45", size = 7),
          axis.title.y.right = element_text(size = 9, angle = 0, vjust = 1))
  ggsave(file.path(cfg$fig_dir, sprintf("pooling_forest_%s.png", gene)), p,
         width = cfg$width, height = 1.6 + 0.27 * nrow(rows), dpi = cfg$dpi)
}
msg("forest plots written to ", cfg$fig_dir)

# ---------------------------------------------------------------------------
# 5. Summary forest plot — all genes x subsets on one page: pooled g, 95% CI
#    (thick) and 95% prediction interval (thin; where the effect in a NEW study
#    is expected to fall). No PI with k = 1; with k = 2 it is extremely wide.
# ---------------------------------------------------------------------------
sf <- pooled[Measure == "g"]
sf[, Subset := factor(Subset, levels = c("All units", "Glomerular", "Tubulointerstitial",
                                         "Whole biopsy"))]
sf[, Gene := factor(Gene, levels = rev(cfg$genes))]
sf[, sig := fifelse(FDR < 0.05, "FDR < 0.05", "n.s.")]
sf[, lab := sprintf("%.2f [%.2f, %.2f]  k=%d%s", estimate, ci_lb, ci_ub, k,
                    fifelse(is.na(I2), "", sprintf(", I²=%.0f%%", I2)))]
p_sf <- ggplot(sf, aes(y = Gene)) +
  geom_vline(xintercept = 0, colour = "grey55", linetype = 2) +
  geom_errorbar(aes(xmin = pi_lb, xmax = pi_ub), width = 0, orientation = "y",
                colour = "grey60", linewidth = 0.4, na.rm = TRUE) +
  geom_errorbar(aes(xmin = ci_lb, xmax = ci_ub), width = 0.25, orientation = "y",
                colour = "#1F3B5C", linewidth = 0.7) +
  geom_point(aes(x = estimate, fill = sig), shape = 23, size = 3, colour = "#1F3B5C") +
  geom_text(aes(x = Inf, label = lab), hjust = 1.02, vjust = -0.9, size = 2.4, colour = "grey30") +
  facet_wrap(~ Subset, nrow = 1) +
  scale_fill_manual(values = c("FDR < 0.05" = "#1B7F79", "n.s." = "white"), name = NULL) +
  coord_cartesian(xlim = c(-4, 4)) +
  labs(title = "Pooled Disease vs Control effect of the RenoRisk panel",
       subtitle = "Hedges' g, DerSimonian–Laird random effects. Thick bar: 95% CI; thin grey bar: 95% prediction interval.",
       x = "Pooled Hedges' g (Disease − Control)", y = NULL,
       caption = paste("BH-FDR across the 5 genes within each column. X axis clipped at ±4;",
                       "Whole biopsy has k = 4 (3 for C4BPA), so its prediction intervals run off the axis.\n",
                       "RNA-Seq units: GSE142025 (whole biopsy), GSE141295 (glomerular), GSE175759 (tubulointerstitial).")) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 8.5),
        plot.caption = element_text(colour = "grey45", size = 7),
        panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        axis.text.y = element_text(face = "bold"), legend.position = "top")
ggsave(file.path(cfg$fig_dir, "pooling_summary_forest.png"), p_sf,
       width = 13, height = 4.2, dpi = cfg$dpi)
msg("wrote ", file.path(cfg$fig_dir, "pooling_summary_forest.png"))

# ---------------------------------------------------------------------------
# 6. Indication-stratified effects (plan §7, Plot 4)
#    Absolute levels are not comparable across cohorts (four are per-gene
#    centred), so the heatmap shows the WITHIN-STUDY effect of each indication
#    against that study's own controls: Hedges' g per batch unit, pooled across
#    units with the same indication by DL random effects (k = 1: the single
#    study estimate). Same min_n rule as the main analysis.
# ---------------------------------------------------------------------------
cfg$indic_file <- file.path(cfg$results_dir, "Stage23", "sample_indication.xlsx")
if (!file.exists(cfg$indic_file))
  stop("missing ", cfg$indic_file, " - run Scripts/stage23_panel_boxplots.R first")
indic <- as.data.table(read_excel(cfg$indic_file, col_types = rep("text", 6)))
# Every cell is pooled: leave out the subsets of GSE30122.
pool_long <- long[!StudyID %in% cfg$subsets]
li <- merge(pool_long, indic[, .(StudyID, GSM, Indication)], by = c("StudyID", "GSM"))
stopifnot(nrow(li) == nrow(pool_long))

ctrl <- li[Indication == "Control", .(m_c = mean(Expression), s_c = sd(Expression), n_c = .N),
           by = .(Gene, Batch, Compartment)]
dis  <- li[Indication != "Control", .(m_d = mean(Expression), s_d = sd(Expression), n_d = .N),
           by = .(Gene, Batch, Compartment, Indication)]
ie <- merge(dis, ctrl, by = c("Gene", "Batch", "Compartment"))[
  n_d >= cfg$min_n & n_c >= cfg$min_n & s_d > 0 & s_c > 0]
es_i <- escalc(measure = "SMD", m1i = m_d, sd1i = s_d, n1i = n_d,
               m2i = m_c, sd2i = s_c, n2i = n_c, data = ie)
ie[, `:=`(g = as.numeric(es_i$yi), g_var = as.numeric(es_i$vi))]

by_ind <- ie[, {
  if (.N >= 2) { f <- rma(yi = g, vi = g_var, data = data.frame(g, g_var), method = cfg$method)
    list(k = .N, n_disease = sum(n_d), estimate = as.numeric(f$beta), ci_lb = f$ci.lb,
         ci_ub = f$ci.ub, p = f$pval, I2 = f$I2,
         units = paste(sort(Batch), collapse = ", ")) }
  else list(k = 1L, n_disease = n_d, estimate = g,
            ci_lb = g - qnorm(0.975) * sqrt(g_var), ci_ub = g + qnorm(0.975) * sqrt(g_var),
            p = 2 * pnorm(-abs(g / sqrt(g_var))), I2 = NA_real_, units = Batch)
}, by = .(Gene, Compartment, Indication)]
by_ind[, FDR := p.adjust(p, "BH")]
write_xlsx(by_ind[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)],
           file.path(cfg$out_dir, "pooling_by_indication.xlsx"))
dropped <- setdiff(unique(li[Indication != "Control", paste(Batch, Indication)]),
                   unique(ie[, paste(Batch, Indication)]))
if (length(dropped))
  msg("indication arms with < ", cfg$min_n, " samples or no in-study controls, left out: ",
      paste(dropped, collapse = "; "))

hm <- copy(by_ind)
hm[, Gene := factor(Gene, levels = cfg$genes)]
hm[, row := sprintf("%s (k=%d, n=%d)", Indication, k, n_disease)]
hm[, row := factor(row, levels = unique(row[order(Compartment, -n_disease)]))]
hm[, row := factor(row, levels = rev(levels(row)))]
hm[, txt := sprintf("%.2f%s", estimate, fifelse(FDR < 0.05, "*", ""))]
lim <- 2
p_hm <- ggplot(hm, aes(Gene, row, fill = pmax(pmin(estimate, lim), -lim))) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = txt), size = 2.7) +
  facet_grid(Compartment ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-lim, lim), name = "Hedges' g\nvs controls") +
  labs(title = "RenoRisk panel by clinical indication",
       subtitle = "Within-study Hedges' g of each indication vs that study's controls, pooled across studies (DL)",
       x = NULL, y = NULL,
       caption = paste0("* BH-FDR < 0.05 across all cells. Colour capped at ±", lim,
                        ". k = studies, n = disease samples. No AKI cohort is in the curated set.")) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 8.5),
        plot.caption = element_text(colour = "grey45", size = 7, hjust = 0),
        axis.text.x = element_text(face = "bold"), panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0))
ggsave(file.path(cfg$fig_dir, "pooling_indication_heatmap.png"), p_hm,
       width = 8.5, height = 1.4 + 0.3 * uniqueN(hm$row), dpi = cfg$dpi)
msg("wrote ", file.path(cfg$fig_dir, "pooling_indication_heatmap.png"))

# ---------------------------------------------------------------------------
# 7. Sensitivity: does the pooled effect depend on the kind of control tissue?
#    Living-donor biopsies are screened healthy kidneys; tumor-nephrectomy
#    tissue is histologically normal but from older patients and taken after
#    surgical ischaemia (relevant to complement genes). Control type per study,
#    from the series-matrix labels / overall design:
#      Living donor      ERCB (GSE104948/GSE104954), GSE108109, GSE200828, GSE116626
#      Tumor nephrectomy GSE99339, GSE129973, GSE96804, GSE121233 (overall design),
#                        GSE142025 (RNA-Seq, whole biopsy: not in the glomerular test),
#                        GSE141295 (RNA-Seq, "normal cortex from renal cell carcinoma
#                        kidney"; glomerular, so in the test)
#      Nephrectomy       GSE175759 ("nephrectomy controls", indication not given;
#                        tubulointerstitial: not in the glomerular test)
#      Not stated        GSE30122 (GEO gives no control source)
#      Deceased donor    GSE111154 (post-mortem organ donors; whole tissue:
#                        not in the glomerular test)
#      Living donor also GSE93798 (glomerular; "Con_L" = living donor, confirmed by
#                        the user 2026-11-04) and GSE35488 (tubulointerstitium:
#                        not in the glomerular test)
#      Reference RNA     GSE66494 ("Normal kidney total RNA", replicates 1-5)
#    ERCB holds both kinds; here its controls are restricted to the living
#    donors, so every unit has ONE control type.
#    Test: glomerular units, Living donor vs Tumor nephrectomy as a moderator
#    (k = 4 vs 5). Glomeruli are the only compartment with enough units of both
#    types, so the test is run there. Not stated and Reference RNA are not
#    tested: each coincides with one disease (GSE30122 family = DKD only,
#    GSE66494 = CKD only), so a difference there would
#    measure disease type.
# ---------------------------------------------------------------------------
CONTROL_TYPE <- c(
  "GSE104948-GPL22945" = "Living donor", "GSE104954-GPL22945" = "Living donor",
  "GSE108109" = "Living donor", "GSE200828" = "Living donor", "GSE116626" = "Living donor",
  "GSE99339-GPL19184" = "Tumor nephrectomy", "GSE129973" = "Tumor nephrectomy",
  "GSE96804" = "Tumor nephrectomy", "GSE121233" = "Tumor nephrectomy",
  "GSE142025" = "Tumor nephrectomy", "GSE141295" = "Tumor nephrectomy",
  "GSE175759" = "Nephrectomy",
  "GSE35488" = "Living donor", "GSE93798" = "Living donor", "GSE111154" = "Deceased donor",
  "GSE30122" = "Not stated",
  "GSE66494" = "Reference RNA")
stopifnot(all(unique(long$StudyID) %in% names(CONTROL_TYPE)))

lc <- merge(long, indic[, .(StudyID, GSM, Label)], by = c("StudyID", "GSM"))
lc[, ControlType := CONTROL_TYPE[StudyID]]
# ERCB: keep only the living donors (blank diagnosis) as controls here.
lc <- lc[!(Condition == "Control" & ControlType == "Living donor" &
           grepl("tumor nephrectom", Label, ignore.case = TRUE))]
cs <- lc[, .(m = mean(Expression), s = sd(Expression), n = .N),
         by = .(Gene, Batch, Compartment, ControlType, Condition)]
cw <- dcast(cs, Gene + Batch + Compartment + ControlType ~ Condition, value.var = c("m", "s", "n"))
cw <- cw[n_Control >= cfg$min_n & n_Disease >= cfg$min_n & s_Control > 0 & s_Disease > 0]
es_c <- escalc(measure = "SMD", m1i = m_Disease, sd1i = s_Disease, n1i = n_Disease,
               m2i = m_Control, sd2i = s_Control, n2i = n_Control, data = cw)
cw[, `:=`(g = as.numeric(es_c$yi), g_var = as.numeric(es_c$vi),
          Pooled = !Batch %in% cfg$subsets)]
write_xlsx(cw[, .(Gene, Batch, Compartment, ControlType, n_Control, n_Disease,
                  g = signif(g, 4), g_var = signif(g_var, 4), Pooled)],
           file.path(cfg$out_dir, "pooling_control_type_study_effects.xlsx"))

ct_pool <- list(); ct_test <- list()
for (gene in cfg$genes) {
  d <- cw[Gene == gene & Pooled]
  # pooled g per control type, glomerular units (k = 1: the study itself)
  for (ct in unique(d[Compartment == "Glomerular", ControlType])) {
    x <- d[Compartment == "Glomerular" & ControlType == ct]
    f <- if (nrow(x) >= 2) rma(yi = g, vi = g_var, data = x, method = cfg$method) else NULL
    ct_pool[[length(ct_pool) + 1]] <- data.table(Gene = gene, ControlType = ct, k = nrow(x),
      estimate = if (is.null(f)) x$g else as.numeric(f$beta),
      ci_lb = if (is.null(f)) x$g - qnorm(0.975) * sqrt(x$g_var) else f$ci.lb,
      ci_ub = if (is.null(f)) x$g + qnorm(0.975) * sqrt(x$g_var) else f$ci.ub,
      I2 = if (is.null(f)) NA_real_ else f$I2)
  }
  # glomerular, Living donor vs Tumor nephrectomy
  ga <- d[Compartment == "Glomerular" & ControlType %in% c("Living donor", "Tumor nephrectomy")]
  ga[, ControlType := factor(ControlType, levels = c("Living donor", "Tumor nephrectomy"))]
  fa <- rma(yi = g, vi = g_var, mods = ~ ControlType, data = ga, method = cfg$method)
  ct_test[[gene]] <- data.table(Gene = gene,
    glom_k_LD = sum(ga$ControlType == "Living donor"),
    glom_k_TN = sum(ga$ControlType == "Tumor nephrectomy"),
    glom_diff_TN_minus_LD = as.numeric(fa$beta[2]), glom_diff_ci_lb = fa$ci.lb[2],
    glom_diff_ci_ub = fa$ci.ub[2], glom_QM_p = fa$QMp)
}
ct_pool <- rbindlist(ct_pool)
ct_test <- rbindlist(ct_test)
ct_test[, glom_FDR := p.adjust(glom_QM_p, "BH")]
write_xlsx(ct_pool[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)],
           file.path(cfg$out_dir, "pooling_control_type_glomerular.xlsx"))
write_xlsx(ct_test[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)],
           file.path(cfg$out_dir, "pooling_control_type_test.xlsx"))
msg("control-type sensitivity (Hedges' g; TN - LD difference in glomerular units):")
print(ct_test[, .(Gene, k_LD = glom_k_LD, k_TN = glom_k_TN,
                  diff = round(glom_diff_TN_minus_LD, 2), p = signif(glom_QM_p, 2),
                  FDR = signif(glom_FDR, 2))])

cp <- copy(ct_pool)
cp[, Gene := factor(Gene, levels = rev(cfg$genes))]
cp[, ControlType := factor(ControlType, levels = c("Living donor", "Tumor nephrectomy", "Not stated"))]
# The test result goes into the gene's axis label, clear of the error bars.
glab <- ct_test[, setNames(sprintf("%s\nTN − LD %.2f, p = %.2f", Gene,
                                   glom_diff_TN_minus_LD, glom_QM_p), Gene)]
p_ct <- ggplot(cp, aes(estimate, Gene, colour = ControlType)) +
  geom_vline(xintercept = 0, colour = "grey55", linetype = 2) +
  geom_errorbar(aes(xmin = ci_lb, xmax = ci_ub), width = 0.2, orientation = "y",
                position = position_dodge(0.6)) +
  geom_point(aes(shape = ControlType), size = 2.6, position = position_dodge(0.6)) +
  scale_y_discrete(labels = glab) +
  scale_colour_manual(values = c("Living donor" = "#1B7F79", "Tumor nephrectomy" = "#B2182B",
                                 "Not stated" = "grey55"), name = "Control tissue") +
  scale_shape_manual(values = c(16, 17, 15), name = "Control tissue") +
  labs(title = "Sensitivity: does the disease effect depend on the control tissue?",
       subtitle = "Glomerular units, pooled Hedges' g (Disease − Control) by control type, DL random effects, 95% CI",
       x = "Pooled Hedges' g", y = NULL,
       caption = paste0("Living donor k = 4 (ERCB, living-donor controls only; GSE108109; GSE200828; GSE93798).  ",
                        "Tumor nephrectomy k = 5 (GSE99339, GSE129973, GSE96804, GSE121233, GSE141295 RNA-Seq).\n",
                        "Not stated k = 1 (GSE30122 glomeruli, DKD only) — shown, not tested.  ",
                        "p = moderator test, Living donor vs Tumor nephrectomy.")) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 8.5),
        plot.caption = element_text(colour = "grey45", size = 7, hjust = 0),
        panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 8.5), legend.position = "top")
ggsave(file.path(cfg$fig_dir, "pooling_control_type_sensitivity.png"), p_ct,
       width = 9, height = 5.5, dpi = cfg$dpi)
msg("wrote ", file.path(cfg$fig_dir, "pooling_control_type_sensitivity.png"))

writeLines(capture.output(sessionInfo()), file.path(cfg$out_dir, "pooling_sessionInfo.txt"))
msg("done")
