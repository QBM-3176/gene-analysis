#!/usr/bin/env Rscript
# =============================================================================
# MetaAnalysis — summary figures that combine the outputs of several steps
#
# SOURCE: CKD_Biomarker_Pooled_Analysis_Plan_Final.pdf — "PROJECT ANALYSIS PLAN &
#   TECHNICAL MASTER PLAN", RenoRisk Validation Across Multiple CKD Indications,
#   Version 2.0 (Final Analysis Plan), July 2026. Project Lead: Amela Jusic;
#   Execution Lead: Stojanko Tudjarski. Quotations below are verbatim.
#
# PLAN REQUIREMENTS TAKEN INTO ACCOUNT HERE
#   Executive summary: "... Stojanko Tudjarski will execute data ingestion,
#     cross-cohort batch harmonization, machine learning modeling, and reporting
#     across 32 public transcriptomic datasets."
#   §4 Cohort Stratification & GEO Dataset Inventory: "The 32 datasets identified
#     in the GEO repository search have been mapped to evaluate RenoRisk target
#     performance across diverse clinical indications: see table in annex 1."
#   Annex 2: "List of the study eligible for meta-analysis." (columns: Condition,
#     GEO Accession, Validated GEO Accession, Cohort / Study, # Samples,
#     # Validated Samples, Sequencing Type, Sample Type, Fits Meta-analysis
#     criteria, Confirmed Quantitative Clinical Variables, GEO Link)
#   §1 Project Roles: "Code generation, pipeline reproducibility, and figure
#     preparation"; "Lead author on executive reports and scientific
#     dissemination".
#   §3 Stage 4 ("Model differential expression of the 5 RenoRisk genes ...") and
#     §5 ("summary-level random-effects meta-analysis as the primary framework")
#     — the two routes compared in the method-agreement figure.
#
# HOW THIS SCRIPT ADDRESSES THEM
#   * Neither figure is a numbered plan plot (§7 Plots 1-5); both support the
#     reporting the plan assigns to the execution lead.
#   * PRISMA-style flow: the selection count comes from the curated dataset
#     inventory (16 GEO series, the same accessions as Annex 2) plus the 8 series
#     added to Data/Raw_filtered on 2026-11-03, not the plan's "32": 24 listed ->
#     7 excluded (2 no control arm, 4 duplicates incl. GSE30528/GSE30529,
#     1 too few patients) -> 17 included (14 microarray + 3 RNA-Seq; GSE141295 and
#     GSE175759 from raw counts via rnaseq_counts_vst.R).
#     GSE137570 / GSE219185 (no control arm) enter the ML inputs only.
#   * Method agreement: the Stage 4 model and the §5 meta-analysis estimate the
#     same Disease - Control contrast by independent routes; their agreement is
#     reported as supporting evidence for both.
#
#   1. Method agreement: Stage 4 (one multi-factor WLS model over all samples,
#      sex-adjusted) vs the pooling strategy (per-study effects pooled by DL
#      random effects). Same quantity — Disease minus Control, log2 — from two
#      independent routes; agreement supports both.
#   2. PRISMA-style flow diagram of dataset selection.
#
# Input : Results/Stage23/panel_uncorrected_long.xlsx, sample_indication.xlsx
#         Results/Stage4/stage4_multifactor_results.xlsx, stage4_per_study_effects.xlsx
#         Results/PoolingStrategy/pooling_estimates.xlsx, pooling_study_effects.xlsx
# Output: Results/Summary/ — tables (.xlsx), method_agreement.png, prisma_flow.png, sessionInfo
#
# Run AFTER stage23_panel_boxplots.R, stage4_MF_stat_mod.R and pooling_strategy.R.
#
# PRISMA counts: 16 datasets from the curated dataset inventory (an Excel sheet
# kept outside this project; the same 16 accessions are in plan Annex 2) + 8
# added to Data/Raw_filtered on 2026-11-03 (by instruction, 3 of them included).
# The plan text's "32 datasets" is not the selection count and is not used.
# Exclusion reasons are curated (cfg$excluded); everything downstream is
# computed from the data.
# =============================================================================

cfg <- list(
  project_dir = path.expand("~/Striki/Quant/quant.2/MetaAnalysis"),
  genes       = c("C4BPA", "COL13A1", "HABP2", "CFH", "CFI"),
  n_listed = 16 + 8,   # curated dataset inventory (Excel) + 8 series added 2026-11-03
# GSE30528/GSE30529 (subsets of GSE30122) were dropped and renamed -dup on 2026-11-03 (by instruction); the guard stays for any future contained dataset.
  subsets  = character(0),   # datasets fully contained in another: not in pooled counts
  rnaseq   = c("GSE142025", "GSE141295", "GSE175759"),   # included RNA-Seq series
  # Exclusion reason per series (see Raw_filtered_inclusion_review.docx)
  excluded = c(GSE137570 = "No control arm (ML inputs only)",
               GSE219185 = "No control arm (ML inputs only)",
               GSE35487  = "Duplicate (same arrays as GSE35488)",
               GSE37460  = "Duplicate (arrays in GSE104948 / GSE99339)",
               GSE30528  = "Duplicate (all samples in GSE30122)",
               GSE30529  = "Duplicate (all samples in GSE30122)",
               GSE1009   = "Too few patients (2 per group)"),
  dpi = 300
)
cfg$results_dir <- file.path(cfg$project_dir, "Results")   # everything the pipeline generates
cfg$out_dir <- file.path(cfg$results_dir, "Summary")
cfg$fig_dir <- cfg$out_dir   # figures sit with the tables

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(readxl); library(writexl)
})
set.seed(1)
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg$fig_dir, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

rd <- function(...) {
  f <- file.path(cfg$results_dir, ...)
  if (!file.exists(f)) stop("missing ", f, " - run the upstream script first")
  as.data.table(read_excel(f))
}

# ---------------------------------------------------------------------------
# 1. Method agreement — Stage 4 vs pooling strategy, log2FC (Disease - Control)
# ---------------------------------------------------------------------------
s4_study <- rd("Stage4", "stage4_per_study_effects.xlsx")
pl_study <- rd("PoolingStrategy", "pooling_study_effects.xlsx")
s4_pool  <- rd("Stage4", "stage4_multifactor_results.xlsx")
pl_pool  <- rd("PoolingStrategy", "pooling_estimates.xlsx")

# Per study unit: Stage 4's sex-adjusted within-unit contrast vs the plain
# difference of means the pooling strategy starts from.
ag_study <- merge(
  s4_study[, .(Gene, Batch, Compartment, s4 = log2FC, s4_lo = log2FC - 1.96 * SE,
               s4_hi = log2FC + 1.96 * SE)],
  pl_study[, .(Gene, Batch, pl = log2FC, pl_lo = log2FC - 1.96 * sqrt(log2FC_var),
               pl_hi = log2FC + 1.96 * sqrt(log2FC_var))],
  by = c("Gene", "Batch"))
ag_study[, Level := "Per study unit"]

# Per compartment: Stage 4's model estimate vs the DL-pooled estimate.
s4c <- s4_pool[Family == "Disease"]
s4c[, Compartment := sub("^Disease vs Control \\| ", "", Term)]
ag_pool <- merge(
  s4c[, .(Gene, Compartment, s4 = log2FC, s4_lo = log2FC - 1.96 * SE, s4_hi = log2FC + 1.96 * SE)],
  pl_pool[Measure == "log2FC" & Subset != "All units",
          .(Gene, Compartment = Subset, pl = estimate, pl_lo = ci_lb, pl_hi = ci_ub)],
  by = c("Gene", "Compartment"))
ag_pool[, Level := "Pooled per compartment"]

ag <- rbind(ag_study, ag_pool, fill = TRUE)
ag[, Gene := factor(Gene, levels = cfg$genes)]
ag[, Level := factor(Level, levels = c("Per study unit", "Pooled per compartment"))]
stats <- ag[, .(n = .N, pearson_r = cor(s4, pl), slope = coef(lm(pl ~ s4))[2],
                same_sign = mean(sign(s4) == sign(pl)), max_abs_diff = max(abs(s4 - pl))),
            by = Level]
write_xlsx(ag, file.path(cfg$out_dir, "method_agreement.xlsx"))
write_xlsx(stats[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)],
           file.path(cfg$out_dir, "method_agreement_stats.xlsx"))
msg("Stage 4 vs pooling strategy agreement:")
print(stats)

facet_lab <- stats[, .(Level, lab = sprintf("%s\nn = %d, r = %.3f, same sign %.0f%%",
                                            Level, n, pearson_r, 100 * same_sign))]
ag[facet_lab, Facet := i.lab, on = "Level"]
ag[, Facet := factor(Facet, levels = facet_lab$lab)]
p1 <- ggplot(ag, aes(s4, pl, colour = Gene, shape = Compartment)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = 2) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_errorbar(data = ag[Level == "Pooled per compartment"],
                aes(ymin = pl_lo, ymax = pl_hi), width = 0, alpha = 0.5) +
  geom_errorbar(data = ag[Level == "Pooled per compartment"],
                aes(xmin = s4_lo, xmax = s4_hi), width = 0, alpha = 0.5, orientation = "y") +
  geom_point(size = 2.2, alpha = 0.85) +
  facet_wrap(~ Facet, nrow = 1, scales = "free") +
  scale_colour_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(title = "Two independent routes to the disease effect agree",
       subtitle = "log2 fold change, Disease − Control. Dashed line: perfect agreement.",
       x = "Stage 4 (multi-factor WLS model, sex-adjusted)",
       y = "Pooling strategy (per-study effects, DL random effects)",
       caption = "Right panel: bars are 95% CIs. Whole biopsy has k = 4 studies (3 for C4BPA), hence the wide pooled intervals.") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 8.5),
        plot.caption = element_text(colour = "grey45", size = 7))
ggsave(file.path(cfg$fig_dir, "method_agreement.png"), p1, width = 11, height = 5.2, dpi = cfg$dpi)
msg("wrote ", file.path(cfg$fig_dir, "method_agreement.png"))

# ---------------------------------------------------------------------------
# 2. PRISMA-style flow diagram
# ---------------------------------------------------------------------------
unc   <- rd("Stage23", "panel_uncorrected_long.xlsx")
indic <- rd("Stage23", "sample_indication.xlsx")
# 14 microarray + 3 RNA-Seq series are included, all in every pooled analysis.
n_inc <- uniqueN(unc$StudyID)
smp   <- unique(unc[!StudyID %in% cfg$subsets, .(StudyID, GSM, Batch, Compartment, Condition)])
comp  <- smp[, .(n = .N, k = uniqueN(Batch)), by = Compartment][order(-n)]
cond  <- smp[, .N, by = Condition]
n_ind <- uniqueN(indic[!StudyID %in% cfg$subsets & Indication != "Control", Indication])
stopifnot(n_inc + length(cfg$excluded) == cfg$n_listed)

ex_n <- table(factor(cfg$excluded, levels = unique(cfg$excluded)))
prisma <- data.table(
  Step  = c("Listed (inventory + added 2026-11-03)", paste("Excluded:", names(ex_n)), "Included"),
  Count = c(cfg$n_listed, as.integer(ex_n), n_inc))
write_xlsx(prisma, file.path(cfg$out_dir, "prisma_counts.xlsx"))

boxes <- data.table(
  x = c(0, 0, 1.3), y = c(3, 1.45, 2.45), w = c(0.95, 0.95, 1.25), h = c(0.55, 1.4, 0.95),
  fill = c("main", "inc", "exc"),
  label = c(
    sprintf("GEO datasets: curated inventory (16)\n+ added 2026-11-03 (8),  n = %d", cfg$n_listed),
    sprintf(paste0("Included in the analysis: n = %d series\n(%d microarray + %d RNA-Seq)\n",
                   "Pooled analyses: %d series, %d units, %d samples\n",
                   "(GSE30122 split by compartment)\n",
                   "%s\n",
                   "Control %d / Disease %d;  %d disease indications"),
            n_inc, n_inc - length(cfg$rnaseq), length(cfg$rnaseq), uniqueN(smp$StudyID), uniqueN(smp$Batch),
            uniqueN(smp$GSM),
            paste(sprintf("%s: %d samples, k = %d", comp$Compartment, comp$n, comp$k), collapse = "\n"),
            cond[Condition == "Control", N], cond[Condition == "Disease", N], n_ind),
    paste0(sprintf("Excluded: n = %d\n", length(cfg$excluded)),
           paste(sprintf("%s: %s", names(ex_n),
                         vapply(names(ex_n), function(r) paste(names(cfg$excluded)[cfg$excluded == r],
                                                               collapse = ", "), "")),
                 collapse = "\n"))))
# Main column: listed -> included. The exclusion box branches off the vertical
# line at its own height.
arrows <- data.table(x = 0, xend = 0, y = 3 - 0.275, yend = 1.45 + 0.7)
side   <- data.table(x = 0, xend = 1.3 - 0.625, y = 2.45, yend = 2.45)
p2 <- ggplot() +
  geom_rect(data = boxes, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2,
                              ymax = y + h / 2, fill = fill), colour = "grey30", linewidth = 0.4) +
  geom_text(data = boxes, aes(x, y, label = label), size = 3, lineheight = 1.05) +
  geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(2.2, "mm"), type = "closed"), linewidth = 0.4) +
  geom_segment(data = side, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(2.2, "mm"), type = "closed"), linewidth = 0.4) +
  scale_fill_manual(values = c(main = "#E8EEF5", inc = "#DCEFEA", exc = "#F6E3E3"),
                    guide = "none") +
  coord_cartesian(xlim = c(-0.55, 2.0), ylim = c(0.65, 3.4), expand = FALSE) +
  labs(title = "Dataset selection (PRISMA-style)") +
  theme_void(base_size = 10) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.background = element_rect(fill = "white", colour = NA))
ggsave(file.path(cfg$fig_dir, "prisma_flow.png"), p2, width = 10.5, height = 4.8, dpi = cfg$dpi)
msg("wrote ", file.path(cfg$fig_dir, "prisma_flow.png"))

writeLines(capture.output(sessionInfo()), file.path(cfg$out_dir, "summary_sessionInfo.txt"))
msg("done")
