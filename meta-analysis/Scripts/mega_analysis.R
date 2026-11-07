#!/usr/bin/env Rscript
# =============================================================================
# MetaAnalysis — mega-analysis within platform tiers (plan §5)
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
#   §3 Stage 3 (study / condition structure the plan models): "... modeling ~
#     StudyID + Condition."
#   §2 Target RenoRisk Biomarker Panel Specification:
#     "The RenoRisk panel consists of 5 specific target biomarkers categorized into
#     primary diagnostic candidates and additional complement system regulators"
#     Main RenoRisk: C4BPA (C4b-Binding Protein Alpha), COL13A1 (Collagen Type XIII
#     Alpha 1 Chain), HABP2 (Hyaluronan Binding Protein 2); Additional RenoRisk:
#     CFH (Complement Factor H), CFI (Complement Factor I).
#
# HOW THIS SCRIPT ADDRESSES THEM (and where it deviates, with the reason)
#   * "mega-analysis": the samples (not summary effects) of every study on a tier
#     are pooled into one model per gene — complementing the summary-level
#     meta-analysis in pooling_strategy.R, which stays the primary framework.
#   * "within platform tiers": tiers are the array manufacturers (Affymetrix,
#     Illumina, Agilent) plus RNA-Seq, so no model mixes platforms. Illumina and
#     Agilent hold one study each; RNA-Seq holds three (GSE142025 whole biopsy,
#     GSE141295 glomerular, GSE175759 tubulointerstitial), one per compartment,
#     so only the Affymetrix tier pools studies within a compartment.
#   * StudyID enters as a random study level and a random disease effect;
#     Condition as a Disease effect per compartment.
#   * The plan does not specify the mega-analysis model; the choices below
#     (mixed model, per-study residual variance, log2 scale) are documented with
#     their reasons.
#
# Plan §5: random-effects meta-analysis is the primary framework, "complemented
# by mega-analysis within platform tiers". The meta-analysis (pooling_strategy.R)
# pools per-study SUMMARY effects; the mega-analysis pools the SAMPLES of all
# studies on one platform tier and fits a single model to them.
#
# Input : Results/Stage23/panel_uncorrected_long.xlsx (written by Stage 2/3)
#         Results/PoolingStrategy/pooling_estimates.xlsx (comparison only)
# Output: Results/MegaAnalysis/ — tables (.xlsx), mega_vs_meta.png, sessionInfo
# Run AFTER stage23_panel_boxplots.R and pooling_strategy.R.
#
# PLATFORM TIERS (from each series' GPL)
#   Affymetrix  12 series, 13 batch units: all 9 glomerular + all 3 tubulointerstitial
#               + GSE111154 (whole tissue, HTA 2.0) (custom-CDF GPL22945/GPL19983/
#               GPL19184/GPL14663, HTA 2.0 GPL17586, U133A2 GPL571)
#   Illumina    GSE116626 (GPL14951), whole biopsy — a single study
#   Agilent     GSE66494  (GPL6480),  whole biopsy — a single study
#   RNA-Seq     GSE142025 (whole biopsy, no C4BPA), GSE141295 (glomerular) and
#               GSE175759 (tubulointerstitial; both DESeq2 VST from counts)
#   Only the Affymetrix tier pools studies within a compartment; elsewhere the
#   mega-analysis IS that study's model (reported with k = 1).
#
# MODEL (per tier x gene; nlme::lme, REML; log2 scale, uncorrected values)
#   Expression ~ 0 + Compartment + Compartment:Disease
#     random  = pdDiag(~ Disease) | Batch   study-specific level and disease
#                                           effect (independent)
#     weights = varIdent(~ 1 | Batch)       residual variance per study (spread
#                                           differs up to ~6.7x between studies)
#   The random intercept absorbs each study's offset, the way removeBatchEffect
#   would, but inside the model. The random disease slope is the sample-level
#   analogue of the meta-analysis tau² (one shared across compartments); if it
#   does not converge the model falls back to a random intercept only (column
#   `model`). One-study tiers: lm(Expression ~ Disease).
#   Disease effects are on the log2 scale, comparable to the meta-analysis
#   log2FC. A standardized (within-study z-score) variant was tried and dropped:
#   z-scoring divides by each study's total SD, which already contains the
#   disease effect, and the between-study disease variance then collapses to 0
#   (falsely precise estimates).
#
# SUBSETS: GSE30528/GSE30529 are subsets of GSE30122. The mega-analysis
# aggregates across datasets, so they are left out (their samples are already
# in GSE30122) and no sample enters twice. This was formerly a sensitivity run
# next to an "all 13 datasets" primary; it is now the only analysis.
# =============================================================================

cfg <- list(
  project_dir = path.expand("~/Striki/Quant/quant.2/MetaAnalysis"),
  genes       = c("C4BPA", "COL13A1", "HABP2", "CFH", "CFI"),
# GSE30528/GSE30529 (subsets of GSE30122) were dropped and renamed -dup on 2026-11-03 (by instruction); the guard stays for any future contained dataset.
  subsets     = character(0),   # datasets fully contained in another: not aggregated
  dpi         = 300
)
cfg$results_dir <- file.path(cfg$project_dir, "Results")   # everything the pipeline generates
cfg$in_file   <- file.path(cfg$results_dir, "Stage23", "panel_uncorrected_long.xlsx")
cfg$meta_file <- file.path(cfg$results_dir, "PoolingStrategy", "pooling_estimates.xlsx")
cfg$out_dir   <- file.path(cfg$results_dir, "MegaAnalysis")
cfg$fig_dir   <- cfg$out_dir   # figures sit with the tables

PLATFORM <- c(
  "GSE104948-GPL22945" = "Affymetrix", "GSE104954-GPL22945" = "Affymetrix",
  "GSE108109" = "Affymetrix", "GSE200828" = "Affymetrix", "GSE99339-GPL19184" = "Affymetrix",
  "GSE121233" = "Affymetrix", "GSE129973" = "Affymetrix", "GSE96804" = "Affymetrix",
  "GSE30122" = "Affymetrix",
  "GSE93798" = "Affymetrix", "GSE35488" = "Affymetrix", "GSE111154" = "Affymetrix",
  "GSE116626" = "Illumina", "GSE66494" = "Agilent", "GSE142025" = "RNA-Seq",
  "GSE141295" = "RNA-Seq", "GSE175759" = "RNA-Seq")

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(nlme); library(readxl); library(writexl)
})
set.seed(1)
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

for (f in c(cfg$in_file, cfg$meta_file))
  if (!file.exists(f)) stop("missing ", f, " - run the upstream script first")

# ---------------------------------------------------------------------------
# 1. Data: uncorrected panel values + platform tier
# ---------------------------------------------------------------------------
# Column types are fixed, not guessed: GSM, StudyID, Condition, Compartment,
# Batch, Gene (text), Expression (numeric).
long <- as.data.table(read_excel(cfg$in_file, col_types = c(rep("text", 6), "numeric")))
long <- long[!StudyID %in% cfg$subsets]   # aggregated analysis: subsets of GSE30122 out
stopifnot(all(unique(long$StudyID) %in% names(PLATFORM)))
long[, Tier := PLATFORM[StudyID]]
long[, Disease := as.integer(Condition == "Disease")]
long[, Compartment := factor(Compartment)]

tiers <- long[, .(units = uniqueN(Batch), studies = uniqueN(StudyID), samples = uniqueN(GSM),
                  compartments = paste(sort(unique(as.character(Compartment))), collapse = ", ")),
              by = Tier]
write_xlsx(tiers, file.path(cfg$out_dir, "mega_tiers.xlsx"))
msg("platform tiers:"); print(tiers)

# ---------------------------------------------------------------------------
# 2. Model per tier x gene
# ---------------------------------------------------------------------------
# A compartment with ONE study in the tier gets that study's own model: in a joint
# lme its random disease slope would be borrowed from the other compartments.
fit_one <- function(d) {
  d <- copy(d); d[, Compartment := droplevels(Compartment)]
  single <- d[, uniqueN(Batch), by = Compartment][V1 == 1, as.character(Compartment)]
  out <- lapply(single, function(cp) {
    f  <- lm(Expression ~ Disease, data = d[Compartment == cp])
    ct <- summary(f)$coefficients
    data.table(Compartment = cp, k = 1L, n = nrow(d[Compartment == cp]),
               estimate = ct["Disease", 1], SE = ct["Disease", 2], df = f$df.residual,
               p = ct["Disease", 4], tau_disease = NA_real_, model = "lm (single study)")
  })
  d <- d[!Compartment %in% single]
  if (!nrow(d)) return(rbindlist(out))
  d[, Compartment := droplevels(Compartment)]
  ctrl <- lmeControl(maxIter = 200, msMaxIter = 200, niterEM = 50, returnObject = FALSE)
  fx   <- Expression ~ 0 + Compartment + Compartment:Disease
  wts  <- varIdent(form = ~ 1 | Batch)
  f <- tryCatch(lme(fx, random = list(Batch = pdDiag(~ Disease)), weights = wts,
                    data = d, method = "REML", control = ctrl),
                error = function(e) NULL)
  model <- "lme, random intercept + disease slope"
  if (is.null(f)) {
    f <- lme(fx, random = ~ 1 | Batch, weights = wts, data = d, method = "REML", control = ctrl)
    model <- "lme, random intercept (slope did not converge)"
  }
  tt  <- summary(f)$tTable
  dis <- grep(":Disease$", rownames(tt))
  cps <- sub("^Compartment(.*):Disease$", "\\1", rownames(tt)[dis])
  vc  <- VarCorr(f)
  rbindlist(c(out, list(data.table(Compartment = cps,
             k = vapply(cps, function(cp) uniqueN(d[Compartment == cp, Batch]), 1L),
             n = vapply(cps, function(cp) nrow(d[Compartment == cp]), 1L),
             estimate = tt[dis, "Value"], SE = tt[dis, "Std.Error"], df = tt[dis, "DF"],
             p = tt[dis, "p-value"],
             tau_disease = if ("Disease" %in% rownames(vc)) as.numeric(vc["Disease", "StdDev"])
                           else NA_real_,
             model = model))))
}

run_mega <- function(dat, analysis) {
  out <- list()
  for (tier in unique(dat$Tier)) for (gene in cfg$genes) if (nrow(dat[Tier == tier & Gene == gene]))
    out[[length(out) + 1]] <- cbind(data.table(Analysis = analysis, Tier = tier, Gene = gene),
                                    fit_one(dat[Tier == tier & Gene == gene]))
  out <- rbindlist(out)
  out[, `:=`(ci_lb = estimate - qt(0.975, df) * SE, ci_ub = estimate + qt(0.975, df) * SE)]
  # BH across the 5 genes, within analysis x tier x compartment.
  out[, FDR := p.adjust(p, "BH"), by = .(Analysis, Tier, Compartment)]
  out[]
}

primary <- sprintf("Primary (%d datasets)", uniqueN(long$StudyID))
mega <- run_mega(long, primary)
setcolorder(mega, c("Analysis", "Tier", "Gene", "Compartment", "k", "n", "estimate", "SE",
                    "ci_lb", "ci_ub", "df", "p", "FDR", "tau_disease", "model"))
num <- names(mega)[vapply(mega, is.numeric, TRUE) & !names(mega) %in% c("k", "n", "df")]
write_xlsx(copy(mega)[, (num) := lapply(.SD, signif, 4), .SDcols = num],
           file.path(cfg$out_dir, "mega_estimates.xlsx"))
msg("mega-analysis, primary (log2, Disease - Control):")
print(mega[Analysis == primary, .(Tier, Gene, Compartment, k, est = round(estimate, 3),
                                  ci_lb = round(ci_lb, 3), ci_ub = round(ci_ub, 3),
                                  FDR = signif(FDR, 2), tau = round(tau_disease, 3))])
if (any(grepl("did not converge", mega$model)))
  msg("random-slope fallback used for: ",
      paste(unique(mega[grepl("did not converge", model), paste(Analysis, Tier, Gene)]),
            collapse = "; "))

# ---------------------------------------------------------------------------
# 3. Mega vs meta: same contrast (log2FC), two pooling routes. In Glomerular and
#    Tubulointerstitial the meta-analysis pools the Affymetrix units plus one
#    RNA-Seq unit each (GSE141295, GSE175759), so the Affymetrix mega is close to,
#    not identical with, the meta's unit set. Whole biopsy: the meta pools Illumina + Agilent +
#    RNA-Seq + Affymetrix GSE111154 (k = 4), the mega reports each tier separately.
# ---------------------------------------------------------------------------
meta <- as.data.table(read_excel(cfg$meta_file))[Measure == "log2FC" & Subset != "All units",
  .(Gene, Compartment = Subset, Tier = "Meta-analysis (all tiers)", k, estimate, ci_lb, ci_ub, FDR)]
cmp <- rbind(mega[Analysis == primary, .(Gene, Compartment, Tier = paste0("Mega: ", Tier), k,
                                          estimate, ci_lb, ci_ub, FDR)],
             meta)
write_xlsx(cmp[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)],
           file.path(cfg$out_dir, "mega_vs_meta.xlsx"))

ag <- merge(mega[Analysis == primary & Tier == "Affymetrix", .(Gene, Compartment, mega = estimate)],
            meta[, .(Gene, Compartment, meta = estimate)], by = c("Gene", "Compartment"))
ag_stats <- ag[, .(n = .N, pearson_r = round(cor(mega, meta), 3),
                   same_sign = mean(sign(mega) == sign(meta)),
                   max_abs_diff = round(max(abs(mega - meta)), 3))]
write_xlsx(ag_stats, file.path(cfg$out_dir, "mega_vs_meta_agreement.xlsx"))
msg("Affymetrix mega vs meta-analysis (Glomerular + Tubulointerstitial):"); print(ag_stats)


# ---------------------------------------------------------------------------
# 4. Figure: mega vs meta per gene x compartment (log2FC)
# ---------------------------------------------------------------------------
pd <- copy(cmp)
pd[, Gene := factor(Gene, levels = rev(cfg$genes))]
pd[, Compartment := factor(Compartment, levels = c("Glomerular", "Tubulointerstitial", "Whole biopsy"))]
pd[, Tier := factor(Tier, levels = c("Meta-analysis (all tiers)", "Mega: Affymetrix",
                                     "Mega: Illumina", "Mega: Agilent", "Mega: RNA-Seq"))]
pd[, sig := fifelse(FDR < 0.05, "FDR < 0.05", "n.s.")]
# group = Tier on every layer: the significance fill would otherwise regroup the
# points and dodge them out of line with their error bars.
p <- ggplot(pd, aes(estimate, Gene, colour = Tier, shape = Tier, group = Tier)) +
  geom_vline(xintercept = 0, colour = "grey55", linetype = 2) +
  geom_errorbar(aes(xmin = ci_lb, xmax = ci_ub), width = 0.25, orientation = "y",
                position = position_dodge(0.6)) +
  geom_point(aes(fill = sig), size = 2.6, position = position_dodge(0.6)) +
  facet_wrap(~ Compartment, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = c("Meta-analysis (all tiers)" = "#1F3B5C", "Mega: Affymetrix" = "#1B7F79",
                                 "Mega: Illumina" = "#B2182B", "Mega: Agilent" = "#E08214",
                                 "Mega: RNA-Seq" = "#7B3294"),
                      name = NULL) +
  scale_shape_manual(values = c(23, 21, 24, 22, 25), name = NULL) +
  scale_fill_manual(values = c("FDR < 0.05" = "grey20", "n.s." = "white"), name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(title = "Mega-analysis within platform tiers vs random-effects meta-analysis",
       subtitle = paste("log2 fold change, Disease − Control, 95% CI. Meta: DL-pooled per-study log2FC;",
                        "mega: one mixed model per tier (random study level and disease effect)"),
       x = "log2 fold change", y = NULL,
       caption = paste0(sprintf("Affymetrix tier: glomerular k = %d, tubulointerstitial k = %d. ",
                                uniqueN(long[PLATFORM[StudyID] == "Affymetrix" & Compartment == "Glomerular", Batch]),
                                uniqueN(long[PLATFORM[StudyID] == "Affymetrix" & Compartment == "Tubulointerstitial", Batch])),
                        "RNA-Seq: GSE141295 (glomerular), GSE175759 (tubulointerstitial), GSE142025 (whole biopsy). ",
                        "Whole biopsy: Affymetrix = GSE111154, Illumina = GSE116626, Agilent = GSE66494.\n",
                        "Filled = FDR < 0.05.")) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 8.5),
        plot.caption = element_text(colour = "grey45", size = 7, hjust = 0),
        panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        axis.text.y = element_text(face = "bold"), legend.position = "top")
ggsave(file.path(cfg$fig_dir, "mega_vs_meta.png"), p, width = 11, height = 4.8, dpi = cfg$dpi)
msg("wrote ", file.path(cfg$fig_dir, "mega_vs_meta.png"))

writeLines(capture.output(sessionInfo()), file.path(cfg$out_dir, "mega_sessionInfo.txt"))
msg("done")
