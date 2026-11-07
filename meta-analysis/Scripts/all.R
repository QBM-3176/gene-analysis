#!/usr/bin/env Rscript
# =============================================================================
# MetaAnalysis — RenoRisk 5-gene panel: per-dataset plots, Stages 2–4 and the
# pooling strategy in one pipeline, plus the ML input workbooks
#
# SOURCE: CKD_Biomarker_Pooled_Analysis_Plan_Final.pdf — "PROJECT ANALYSIS PLAN &
#   TECHNICAL MASTER PLAN", RenoRisk Validation Across Multiple CKD Indications,
#   Version 2.0 (Final Analysis Plan), July 2026. Quotations are verbatim.
#
# PLAN REQUIREMENTS AND WHERE THEY ARE IMPLEMENTED
#   §3 Stage 2 — Target Gene Mapping (Part B):
#     "Specifically extract probe sets/transcripts corresponding to C4BPB,
#     COL13A1, HABP2, CFH, and CFI. For microarrays with multiple probes per gene,
#     collapse using WGCNA::collapseRows (selecting maximum mean expression)."
#   §3 Stage 3 — Empirical Bayes Batch Correction (Part D):
#     "Apply ComBat-seq for raw integer RNA-Seq counts and limma::removeBatchEffect
#     / sva::ComBat for normalized continuous expression, modeling ~ StudyID +
#     Condition."
#   §3 Stage 4 — Multi-Factor Statistical Modeling (Part E):
#     "Model differential expression of the 5 RenoRisk genes controlling for
#     tissue compartment (glomerular vs. tubulointerstitial) and patient clinical
#     covariates."
#   §5 Methodology: Pooling Strategy (Part F):
#     "Per-Study Effect Size: ... compute Hedges' g standardized mean difference
#     (or log2 fold change) and variance in each dataset." "Pooled Summary
#     Estimate: Combine study effect sizes using DerSimonian-Laird random-effects
#     models via the metafor package in R." "Heterogeneity Evaluation: Evaluate
#     Cochran's Q and I² statistics ..."
#   §5 "... complemented by mega-analysis within platform tiers." (Part G)
#   §4 Cohort Stratification & GEO Dataset Inventory — PRISMA flow (Part H).
#   §7 Plot 1 (PCA / UMAP pre- vs post-correction) — Part D.
#   §7 Plot 2 (forest plots for the 5 genes)         — Part F.
#   §7 Plot 4 (indication-stratified heatmap)        — Part F.
#   Per-dataset boxplots and heatmaps (no plan section; one figure set per
#   dataset, drawn standalone)                       — Part C.
#   Stage 4 vs pooling method agreement (reporting)  — Part H.
#
# PIPELINE
#   Part A  setup: configuration, packages, the one cohort registry
#   Part B  load every series matrix ONCE; map probes to genes; MaxMean collapse
#           (Stage 2). Panel genes and the sex-marker genes share this layer.
#   Part C  per-dataset boxplots (Wilcoxon vs reference arm) and heatmaps
#   Part D  Stage 3: one limma::removeBatchEffect model over all samples;
#           harmonised figures; PCA/UMAP before vs after; residual spread
#   Part E  Stage 4: sex (reported or inferred from XIST / Y genes); one WLS
#           model per gene, Expression ~ 0 + Batch + Disease:Compartment + Sex
#   Part F  pooling strategy: Hedges' g / log2FC per batch unit, DL random
#           effects (overall, per compartment, moderator), leave-one-out,
#           forest plots, indication heatmap, control-tissue sensitivity
#   Part G  mega-analysis: one nlme mixed model per platform tier x gene
#           (random study level + disease effect), comparison with the
#           meta-analysis
#   Part H  summary: Stage 4 vs pooling agreement; PRISMA-style selection flow
#   Part I  ML input: one row per GSM (1322, SUBSETS left out; incl. the two
#           RNA-Seq series without controls, GSE137570-GPL23227 and GSE219185,
#           which enter the ML inputs only), five genes, Sex, Platform,
#           Disease_Label; all.xlsx (uncorrected) and all_bn.xlsx (corrected)
#   Parts D–I pass data in memory; every hand-off table is still written to
#   Results/ as .xlsx, where it is a deliverable.
#
# DEVIATIONS FROM THE PLAN (reasons)
#   * Stage 1 (raw re-quantification) is retired: inputs are the submitter-
#     processed series-matrix values (five pipelines). RNA-Seq joins the same
#     bag (by instruction, 2026-11-03): GSE142025 is a regular cohort, its values
#     the submitter's normalised log2 table (<id>_expression.txt, gene symbols;
#     C4BPA absent), so ComBat-seq never applies. GSE141295 and GSE175759 are
#     deposited as counts: DESeq2 VST (rnaseq_counts_vst.R), plots via DESeq2.
#   * Panel gene C4BPB -> C4BPA: replaced at the request of the Project Lead.
#   * Stage 3 follows the plan, "~ StudyID + Condition" (by instruction,
#     2026-11-04): batch = StudyID, design = intercept + Disease; ComBat-seq on
#     the two count-deposited RNA-Seq series before their VST. Compartment is
#     aliased with study (except GSE30122), so those level differences are removed.
#   * Stage 4: no cohort carries eGFR, %TIF, proteinuria or age; sex is the only
#     patient covariate. Weights = 1 / residual variance per batch unit.
#   * Pooling: unit = batch unit; Hedges' g primary (scale-free: four cohorts are
#     per-gene centred). Plot 4 has no AKI or CKD-stage rows (no such data).
#
# STANDING CAVEATS
#   * GSE30528 and GSE30529 (all 44 samples in GSE30122) were dropped on
#     2026-11-03 (by instruction; folders renamed -dup). A dataset fully contained
#     in another (`contained_in`, SUBSETS — now empty) would be used in single-
#     dataset calculations only. Aggregated: 17 series, 18 units, 1078 samples.
#   * Four cohorts are per-gene centred and one is linear (log2-transformed
#     here); GSE142025 is RNA-Seq (submitter-normalised log2), GSE141295 and
#     GSE175759 RNA-Seq DESeq2 VST made from their deposited counts by
#     Scripts/by-claude/rnaseq_counts_vst.R (run it first). Absolute expression is not
#     comparable across all cohorts.
#   * A panel gene missing from a deposit (C4BPA in GSE142025) stays in the
#     panel: that cohort's values are NA and every per-gene step uses the
#     samples that have the gene (Stage 3 is fitted gene by gene).
#
# Input : Data/Raw_filtered/<GSE>/<id>_series_matrix.txt (14 series matrices);
#         <id>_expression.txt for the RNA-Seq series (GSE142025; Part I also
#         GSE137570-GPL23227 and GSE219185)
# Output: Results/SingleDatasetPlots/, Stage23/, Stage4/, PoolingStrategy/,
#         MegaAnalysis/, Summary/, ML_training/ — tables (.xlsx), figures,
#         <name>_sessionInfo.txt
# Usage : Rscript Scripts/all.R   (no arguments, no network; ~1 min)
# =============================================================================


# #############################################################################
# PART A — SETUP
# #############################################################################

cfg <- list(
  project_dir = path.expand("~/Striki/Quant/quant.2/MetaAnalysis"),
  genes       = c("C4BPA", "COL13A1", "HABP2", "CFH", "CFI"),
  dpi         = 300,
  plots = list(                # Part C — per-dataset plots
    min_n       = 2,           # a group needs this many samples to enter a Wilcoxon
    heatmap_top = 50,          # most variable probes per dataset heatmap
    width = 8, height = 6),
  s23 = list(                  # Part D — Stage 3
    min_group_n = 3,           # drop a batch arm with fewer samples than this
    gridlines   = "vertical",  # "vertical" | "horizontal" | "both" | "none"
    width = 11, height = 4.875,          # by-compartment figure
    pc_width = 13, pc_height = 22,       # per-cohort figure
    pca_width = 12, pca_height = 5.5,    # PCA and UMAP figures
    umap_neighbors = 15, umap_min_dist = 0.1),
  s4 = list(                   # Part E — Stage 4
    x_gene      = "XIST",
    y_genes     = c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY"),
    min_sex_sep = 4,           # inferred sex used only if the two clusters are this
                               # many within-cluster SDs apart; else left unknown
    width = 9, height = 5),
  pool = list(                 # Part F — pooling strategy
    min_n  = 3,                # per arm, for a unit to enter a meta-analysis
    method = "DL",             # DerSimonian–Laird, per the plan
    width  = 9),
  mega = list(                 # Part G — mega-analysis
    width = 11, height = 4.8),
  summary = list(              # Part H — summary figures
    n_listed = 16 + 8,         # curated inventory (Excel, = plan Annex 2) + 8 added 2026-11-03
    # Exclusion reason per series (see Raw_filtered_inclusion_review.docx)
    excluded = c(GSE137570 = "No control arm (ML inputs only)",
                 GSE219185 = "No control arm (ML inputs only)",
                 GSE35487  = "Duplicate (same arrays as GSE35488)",
                 GSE37460  = "Duplicate (arrays in GSE104948 / GSE99339)",
                 GSE30528  = "Duplicate (all samples in GSE30122)",
                 GSE30529  = "Duplicate (all samples in GSE30122)",
                 GSE1009   = "Too few patients (2 per group)"))
)
cfg$raw_dir     <- file.path(cfg$project_dir, "Data", "Raw_filtered")
cfg$results_dir <- file.path(cfg$project_dir, "Results")
cfg$out <- list(plots = file.path(cfg$results_dir, "SingleDatasetPlots"),
                s23   = file.path(cfg$results_dir, "Stage23"),
                s4    = file.path(cfg$results_dir, "Stage4"),
                pool  = file.path(cfg$results_dir, "PoolingStrategy"),
                mega  = file.path(cfg$results_dir, "MegaAnalysis"),
                summary = file.path(cfg$results_dir, "Summary"),
                ml    = file.path(cfg$results_dir, "ML_training"),
                rnaseq = file.path(cfg$results_dir, "RNASeq"),   # as rnaseq_panel.R
                ml_input = file.path(cfg$project_dir, "Data", "ML_input"))   # all.xlsx
for (d in cfg$out) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Attached: the packages used unqualified. Everything else is called with ::
# (WGCNA, limma, AnnotationDbi, ggsignif, pheatmap, uwot, cluster, nlme), so no
# package masks a base function (WGCNA would mask cor()).
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(metafor); library(writexl)
})
msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)
write_session <- function(dir, name)
  writeLines(capture.output(sessionInfo()), file.path(dir, paste0(name, "_sessionInfo.txt")))
t_start <- Sys.time()

# Entrez IDs: the custom-CDF cohorts address probes as "<entrez>_at".
# Panel genes first, then the sex markers used in Stage 4.
ENTREZ <- c(C4BPA = "722", COL13A1 = "1305", HABP2 = "3026", CFH = "3075", CFI = "3426",
            XIST = "7503", RPS4Y1 = "6192", DDX3Y = "8653", KDM5D = "8284",
            UTY = "7404", EIF1AY = "9086")
SEX_GENES <- c(cfg$s4$x_gene, cfg$s4$y_genes)

# Array vendor per probe-ID route: shown in the per-dataset plot subtitles, and
# the platform tier of the mega-analysis. Not every microarray is Affymetrix
# (GSE116626 Illumina, GSE66494 Agilent).
VENDOR <- c(entrez = "Affymetrix", hta20transcriptcluster.db = "Affymetrix",
            hgu133a2.db = "Affymetrix", illuminaHumanv4.db = "Illumina",
            hgug4112a.db = "Agilent", symbol = "RNA-Seq")

# -----------------------------------------------------------------------------
# The cohort registry — the one place every per-cohort fact lives.
#   id            series-matrix basename (<id>_series_matrix.txt)
#   dir           Data/Raw_filtered subfolder (also the per-dataset output folder)
#   idtype        probe-ID route: "entrez" (custom CDF), a ChipDb package, or
#                 "symbol" (RNA-Seq: rows are gene symbols already)
#   expr_file     RNA-Seq: the series-matrix table is empty; values come from
#                 this file in the dataset folder (header Symbol, GSM columns)
#   compartment   tissue, or "split" (GSE30122 holds both, read per sample)
#   cond_field    where the condition label lives: "source", "description" or a
#                 characteristics key; `control` is a regex on it -> Control
#   blank_control ERCB: an EMPTY diagnosis is a living donor (Control); its
#                 tumor nephrectomies are Control too (by instruction, 2026-10-09)
#   log2          deposited on a linear scale (Agilent) -> log2(x + 1)
#   plot_field    grouping field for the per-dataset plots (default cond_field)
#   blank_label   label for an empty grouping value in the per-dataset plots
#   ref           reference arm of the per-dataset significance brackets
#   y_lab         y-axis label of the per-dataset plots (scale as deposited)
#   control_type  control tissue, for the pooling sensitivity check
#   title_sex     sex is recorded only in the sample title ("..._M_67_..."): reported sex
#   contained_in  the dataset holding all of this one's samples: the dataset is
#                 used in single-dataset calculations, never in aggregated ones
# Order = the order samples enter the pooled tables.
# -----------------------------------------------------------------------------
COHORTS <- list(
  list(id = "GSE104948-GPL22945", dir = "GSE104948", idtype = "entrez",
       compartment = "Glomerular", cond_field = "diagnosis", blank_control = TRUE,
       control = "tumor nephrectom", blank_label = "LD", ref = "LD",
       y_lab = "RMA-normalised expression (log2)", control_type = "Living donor"),
  list(id = "GSE104954-GPL22945", dir = "GSE104954", idtype = "entrez",
       compartment = "Tubulointerstitial", cond_field = "diagnosis", blank_control = TRUE,
       control = "tumor nephrectom", blank_label = "LD", ref = "LD",
       y_lab = "RMA-normalised expression (log2)", control_type = "Living donor"),
  list(id = "GSE108109", dir = "GSE108109", idtype = "entrez",
       compartment = "Glomerular", cond_field = "source", control = "living donor",
       ref = "LD", y_lab = "Normalised expression (log2)", control_type = "Living donor"),
  list(id = "GSE200828", dir = "GSE200828", idtype = "entrez",
       compartment = "Glomerular", cond_field = "source", control = "living donor",
       ref = "LD", y_lab = "RMA-normalised expression (log2)", control_type = "Living donor"),
  list(id = "GSE99339-GPL19184", dir = "GSE99339", idtype = "entrez",
       compartment = "Glomerular", cond_field = "description", control = "tumor nephrectomy",
       ref = "TN", y_lab = "Normalised expression (log2)", control_type = "Tumor nephrectomy"),
  # Added 2026-11-03 (by instruction): GSE93798 (IgAN glomeruli, custom CDF v19),
  # GSE35488 (IgAN tubulointerstitium, Entrez CDF v10: bare Entrez IDs),
  # GSE111154 (post-mortem cortex, early DN vs non-diabetic donors, HTA 2.0).
  list(id = "GSE93798", dir = "GSE93798", idtype = "entrez",
       compartment = "Glomerular", cond_field = "source", control = "control",
       ref = "Control", y_lab = "RMA-normalised expression (log2)",
       control_type = "Living donor", title_sex = TRUE),   # "Con_L" = living donor (user, 2026-11-04)
  list(id = "GSE35488", dir = "GSE35488", idtype = "entrez",
       compartment = "Tubulointerstitial", cond_field = "source", control = "living donor",
       ref = "LD", y_lab = "RMA-normalised expression (log2)", control_type = "Living donor"),
  list(id = "GSE111154", dir = "GSE111154", idtype = "hta20transcriptcluster.db",
       compartment = "Whole biopsy", cond_field = "source", control = "control",
       ref = "Control", y_lab = "Normalised expression (log2, Partek)",
       control_type = "Deceased donor"),   # post-mortem organ donors
  list(id = "GSE121233", dir = "GSE121233", idtype = "hta20transcriptcluster.db",
       compartment = "Glomerular", cond_field = "disease state", control = "^control",
       ref = "Control", y_lab = "RMA-normalised expression (log2)",
       control_type = "Tumor nephrectomy"),   # per its overall design
  list(id = "GSE129973", dir = "GSE129973", idtype = "hta20transcriptcluster.db",
       compartment = "Glomerular", cond_field = "source", control = "tumor nephrectom",
       ref = "TN", y_lab = "RMA-normalised expression (log2)", control_type = "Tumor nephrectomy"),
  list(id = "GSE96804", dir = "GSE96804", idtype = "hta20transcriptcluster.db",
       compartment = "Glomerular", cond_field = "source", control = "tumor nephrectom",
       ref = "TN", y_lab = "RMA-normalised expression (log2)", control_type = "Tumor nephrectomy"),
  # GSE30122 carries both compartments; the per-dataset brackets compare every
  # group with Ctrl-Glom, as the other datasets compare with one reference.
  list(id = "GSE30122", dir = "GSE30122", idtype = "hgu133a2.db",
       compartment = "split", cond_field = "source", control = "control kidney",
       ref = "Ctrl-Glom", y_lab = "Per-gene centred expression (log2 ratio)",
       control_type = "Not stated"),
  list(id = "GSE116626", dir = "GSE116626", idtype = "illuminaHumanv4.db",
       compartment = "Whole biopsy", cond_field = "subject status", control = "living donor",
       plot_field = "source", ref = "LD", y_lab = "Per-gene centred expression (log2 ratio)",
       control_type = "Living donor"),
  list(id = "GSE66494", dir = "GSE66494", idtype = "hgug4112a.db", log2 = TRUE,
       compartment = "Whole biopsy", cond_field = "sample type", control = "normal kidney",
       ref = "Normal", y_lab = "log2 expression (Agilent, linear source)",
       control_type = "Reference RNA"),
  # RNA-Seq (by instruction, 2026-11-03: same bag as the microarrays). The series
  # matrix holds labels only; values are the submitter's normalised log2 table in
  # `expr_file`, rows addressed by gene symbol. C4BPA is not in this deposit.
  list(id = "GSE142025", dir = "GSE142025", idtype = "symbol",
       expr_file = "GSE142025_expression.txt",
       compartment = "Whole biopsy", cond_field = "group", control = "^control$",
       ref = "Control", y_lab = "Submitter-normalised log2 (RNA-Seq)",
       control_type = "Tumor nephrectomy"),
  # RNA-Seq deposited as raw counts only (added 2026-11-03, by instruction): values
  # are DESeq2 VST made by Scripts/by-claude/rnaseq_counts_vst.R from GEO <id>_RAW.tar (by
  # instruction, 2026-11-04); `counts = TRUE` sends their plots to the DESeq2 branch.
  list(id = "GSE141295", dir = "GSE141295", idtype = "symbol",
       expr_file = "GSE141295_expression.txt", counts = TRUE,
       compartment = "Glomerular", cond_field = "tissue1", control = "normal cortex",
       ref = "TN", y_lab = "DESeq2 VST (RNA-Seq, from counts)",
       control_type = "Tumor nephrectomy"),   # normal cortex of renal-cell-carcinoma kidneys
  list(id = "GSE175759", dir = "GSE175759", idtype = "symbol",
       expr_file = "GSE175759_expression.txt", counts = TRUE,
       compartment = "Tubulointerstitial", cond_field = "diagnosis", control = "^control$",
       ref = "Control", y_lab = "DESeq2 VST (RNA-Seq, from counts)",
       control_type = "Nephrectomy")          # "nephrectomy controls", indication not given
  # GSE30528/GSE30529 (all samples in GSE30122) were dropped and their folders renamed
  # -dup on 2026-11-03 (by instruction). `contained_in` stays for any future case.
)
names(COHORTS) <- vapply(COHORTS, `[[`, "", "id")
CONTROL_TYPE <- vapply(COHORTS, `[[`, "", "control_type")
# Datasets left out of every aggregated calculation (contained in another one).
SUBSETS <- names(Filter(function(co) !is.null(co$contained_in), COHORTS))
stopifnot(all(vapply(COHORTS[SUBSETS], `[[`, "", "contained_in") %in% names(COHORTS)))

# Specific diagnosis ("indication") per sample, for plan Plot 4. Matched on the
# condition field; first match wins, so mixed and negated labels come first.
# Controls are "Control" whatever their wording; an unmatched Disease label is
# an error, so a new cohort cannot silently fall into the wrong row.
INDICATIONS <- c(
  "focal segmental glomerul[a-z]* sclerosis/minimal|glomerulosclerosis and minimal" = "FSGS/MCD",
  "diabetic|dkd|_dn$"         = "DKD",
  "focal|fsgs"                = "FSGS",
  "minimal change"            = "MCD",
  "membranous"                = "MN",
  "anca"                      = "ANCA vasculitis",
  "lupus"                     = "Lupus nephritis",
  "non-igan"                  = "Other GN",
  "iga|immunoglobulin a nephropathy" = "IgAN",
  "hypertensive"              = "Hypertensive nephropathy",
  "thin membrane"             = "Thin membrane disease",
  "glomerular disease, other" = "Other GN",
  "tumor nephrectom"          = "Tumor nephrectomy",
  "renal biopsy specimens"    = "CKD (unspecified)"
)

# Short group labels for the per-dataset plots — one vocabulary for all
# cohorts. Order matters: the FSGS+MCD compound before plain FSGS.
ABBREV <- list(
  c("glomerulus of control kidney",                           "Ctrl-Glom"),
  c("tubuli of control kidney",                               "Ctrl-Tub"),
  c("glomerulus of dkd kidney",                               "DKD-Glom"),
  c("tubuli of dkd kidney",                                   "DKD-Tub"),
  c("focal segmental.*(minimal change|mcd)",                  "FSGS+MCD"),
  c("focal (segmental|and segmental).*(sclerosis|scelrosis)", "FSGS"),
  c("minimal_lesion_igan",                                    "IgAN-min"),
  c("chronic_lesion_igan",                                    "IgAN-chr"),
  c("mixed_lesions_igan",                                     "IgAN-mix"),
  c("active_lesion_igan",                                     "IgAN-act"),
  c("non.?igan gn",                                           "non-IgAN"),
  c("iga nephropathy|from human igan$|immunoglobulin a nephropathy", "IgAN"),
  c("minimal change disease",                                 "MCD"),
  c("anca.?associated vasculitis",                            "RPGN"),
  c("membranous (nephropathy|glomerulonephritis)",            "MN"),
  c("systemic lupus|lupus nephritis",                          "SLE"),
  c("hypertensive nephropathy",                               "HTN"),
  c("thin membrane disease",                                  "TMD"),
  c("early diabetic nephropathy",                             "Early_DN"),
  c("non-diabetic control|from human control$",               "Control"),
  c("diabetic nephropathy",                                   "DN"),
  c("tumor nephrectom|normal cortex from renal cell carcinoma", "TN"),
  c("living donor|^ld$",                                      "LD"),
  c("glomerular disease, other",                              "Other"),
  c("normal kidney",                                          "Normal"),
  c("renal biopsy specimens",                                 "CKD"),
  c("^control$",                                              "Control")
)
GROUP_PALETTE <- c(
  LD = "#4DAF4A", Control = "#4DAF4A", Normal = "#4DAF4A", TN = "#66C2A5",
  `Ctrl-Glom` = "#4DAF4A", `Ctrl-Tub` = "#1B9E77",
  DN = "#E41A1C", `DKD-Glom` = "#E41A1C", `DKD-Tub` = "#B2182B", CKD = "#E41A1C",
  FSGS = "#377EB8", `FSGS+MCD` = "#984EA3", MCD = "#FF7F00", RPGN = "#A65628",
  MN = "#F781BF", SLE = "#999999", HTN = "#8DA0CB", TMD = "#E7298A",
  IgAN = "#377EB8", `IgAN-min` = "#80B1D3", `IgAN-chr` = "#FB8072",
  `IgAN-mix` = "#BEBADA", `IgAN-act` = "#FDB462", `non-IgAN` = "#B3DE69",
  Other = "#CCCCCC", Early_DN = "#FB8072", Advanced_DN = "#E41A1C"
)
COND_PALETTE <- c(Control = "#4DAF4A", Disease = "#E41A1C")


# #############################################################################
# PART B — LOAD AND MAP (Stage 2: target gene mapping)
# #############################################################################

# Series matrix -> sample metadata + expression table. No GPL fetch. An empty
# table (RNA-Seq series) gives a 0 x 0 matrix; the loader then reads `expr_file`.
read_series_matrix <- function(path) {
  lines <- readLines(path, warn = FALSE)
  meta  <- grep("^!Sample_", lines, value = TRUE)
  parts <- lapply(meta, function(ln) {
    p <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    list(key = sub("^!Sample_", "", p[1]), vals = gsub('^"|"$', "", p[-1]))
  })
  keys  <- make.unique(vapply(parts, `[[`, "", "key"))
  pheno <- as.data.frame(setNames(lapply(parts, `[[`, "vals"), keys), stringsAsFactors = FALSE)
  rownames(pheno) <- pheno$geo_accession
  b <- which(trimws(lines) == "!series_matrix_table_begin") + 1
  e <- which(trimws(lines) == "!series_matrix_table_end")   - 1
  expr <- if (e >= b)
    as.matrix(read.table(text = paste(lines[b:e], collapse = "\n"), header = TRUE,
                         sep = "\t", row.names = 1, check.names = FALSE, quote = '"',
                         comment.char = ""))
  else matrix(numeric(0), 0, 0)
  list(pheno = pheno, expr = expr)
}

# Raw values of a metadata field: "source", "description", or a
# characteristics key (the "key: " prefix removed). Case is preserved.
field_values <- function(pheno, field) {
  if (field == "source")      return(pheno$source_name_ch1)
  if (field == "description") return(pheno$description)
  cols <- grep("^characteristics_ch1", names(pheno), value = TRUE)
  hit  <- cols[vapply(cols, function(c)
            any(grepl(paste0("^", field, ":"), pheno[[c]], ignore.case = TRUE)), TRUE)]
  if (!length(hit)) stop("field not found: ", field)
  sub(paste0("^", field, ":\\s*"), "", pheno[[hit[1]]], ignore.case = TRUE)
}

derive_condition <- function(pheno, co) {
  v <- tolower(field_values(pheno, co$cond_field))
  ctrl <- grepl(co$control, v, ignore.case = TRUE)
  if (isTRUE(co$blank_control)) ctrl <- ctrl | trimws(v) == ""
  factor(ifelse(ctrl, "Control", "Disease"), levels = c("Control", "Disease"))
}

derive_indication <- function(pheno, co, condition) {
  v <- tolower(field_values(pheno, co$cond_field))
  out <- rep(NA_character_, length(v))
  for (pat in names(INDICATIONS))
    out[is.na(out) & grepl(pat, v, ignore.case = TRUE)] <- INDICATIONS[[pat]]
  out[condition == "Control"] <- "Control"
  if (anyNA(out)) stop(co$id, ": no indication for label(s): ",
                       paste(unique(v[is.na(out)]), collapse = " | "))
  data.table(Label = trimws(v), Indication = out)
}

derive_compartment <- function(pheno, co) {
  if (co$compartment != "split") return(rep(co$compartment, nrow(pheno)))
  s <- tolower(pheno$source_name_ch1)       # GSE30122 interleaves both
  ifelse(grepl("glomerul", s), "Glomerular",
         ifelse(grepl("tubul", s), "Tubulointerstitial", NA_character_))
}

# Per-dataset plot groups: abbreviated labels from the grouping field.
group_labels <- function(pheno, co) {
  raw <- trimws(field_values(pheno, if (is.null(co$plot_field)) co$cond_field else co$plot_field))
  if (!is.null(co$blank_label)) raw[raw == ""] <- co$blank_label
  v <- tolower(trimws(raw)); out <- rep(NA_character_, length(v))
  for (rule in ABBREV) out[is.na(out) & grepl(rule[1], v)] <- rule[2]
  ifelse(is.na(out), raw, out)             # unmatched labels pass through verbatim
}

# Probe -> SYMBOL for the whole array, once per non-Entrez cohort. The ChipDb
# object is exported under the full package name, ".db" suffix included.
probe_symbols_all <- function(idtype, probes) {
  if (idtype %in% c("entrez", "symbol")) return(NULL)
  if (!requireNamespace(idtype, quietly = TRUE))
    stop("annotation package not installed: ", idtype)
  suppressMessages(AnnotationDbi::mapIds(getExportedValue(idtype, idtype), keys = probes,
                                         column = "SYMBOL", keytype = "PROBEID",
                                         multiVals = "first"))
}

# Probes carrying `genes`: a named vector probe -> gene symbol.
gene_probes <- function(ck, genes) {
  probes <- rownames(ck$sm$expr)
  if (ck$co$idtype == "symbol") {   # RNA-Seq deposits: rows are gene symbols already
    hit <- intersect(genes, probes)
    return(setNames(hit, hit))
  }
  if (ck$co$idtype == "entrez") {   # "<entrez>_at", or bare "<entrez>" (GSE35488, CDF v10)
    want <- c(paste0(ENTREZ[genes], "_at"), ENTREZ[genes])
    hit  <- intersect(want, probes)
    return(setNames(rep(genes, 2)[match(hit, want)], hit))
  }
  sym <- ck$sym[!is.na(ck$sym) & ck$sym %in% genes]
  setNames(unname(sym), names(sym))
}

# Multi-probe genes -> the probe with the highest mean (plan: collapseRows,
# MaxMean). Rows of the result are named by gene symbol.
collapse_genes <- function(x, map)
  WGCNA::collapseRows(x[names(map), , drop = FALSE], rowGroup = unname(map),
                      rowID = names(map), method = "MaxMean")$datETcollapsed

on_log2 <- function(x, co) if (isTRUE(co$log2)) log2(x + 1) else x   # Agilent linear -> log2

msg("Part B — loading ", length(COHORTS), " series matrices")
CK <- list()   # per cohort: registry entry, series matrix, symbol map, panel matrix
for (co in COHORTS) {
  f <- file.path(cfg$raw_dir, co$dir, paste0(co$id, "_series_matrix.txt"))
  if (!file.exists(f)) { warning(co$id, ": series matrix not found"); next }
  ck <- list(co = co, sm = read_series_matrix(f))
  if (!is.null(co$expr_file)) {   # RNA-Seq: the series-matrix table is empty
    ex <- fread(file.path(cfg$raw_dir, co$dir, co$expr_file))
    ck$sm$expr <- as.matrix(ex[, -1]); rownames(ck$sm$expr) <- ex[[1]]
    # The file may omit samples (GSE175759: submitter-flagged technical outliers);
    # series-matrix samples without values are dropped.
    stopifnot(all(colnames(ck$sm$expr) %in% rownames(ck$sm$pheno)))
    ck$sm$pheno <- ck$sm$pheno[rownames(ck$sm$pheno) %in% colnames(ck$sm$expr), , drop = FALSE]
    ck$sm$expr <- ck$sm$expr[, rownames(ck$sm$pheno)]
  }
  ck$rnaseq <- any(grepl("RNA-Seq", unlist(ck$sm$pheno[, grep("^library_strategy",
                         names(ck$sm$pheno)), drop = FALSE]), ignore.case = TRUE))
  ck$sym <- tryCatch(probe_symbols_all(co$idtype, rownames(ck$sm$expr)),
                     error = function(e) { warning(co$id, ": ", conditionMessage(e)); NULL })
  map <- if (nrow(ck$sm$expr)) gene_probes(ck, cfg$genes) else character(0)
  if (length(map)) {
    x <- ck$sm$expr[names(map), , drop = FALSE]
    if (isTRUE(co$log2)) stopifnot(min(x, na.rm = TRUE) > -1)
    ck$panel <- collapse_genes(on_log2(x, co), map)
  }
  CK[[co$id]] <- ck
  msg(co$id, ": ", nrow(ck$sm$expr), " probes, ", nrow(ck$sm$pheno), " samples, ",
      if (is.null(ck$panel)) "no panel genes" else paste(nrow(ck$panel), "panel genes"))
}


# #############################################################################
# PART C — PER-DATASET BOXPLOTS AND HEATMAPS
#   Each dataset drawn standalone on its deposited scale: one boxplot per gene
#   (groups from the dataset's own labels, Wilcoxon vs the reference arm) and
#   two heatmaps (top-N most variable probes; the panel genes found).
#   GSE142025 (RNA-Seq, normalised log2 from `expr_file`) takes the microarray
#   path. The counts branch (DESeq2 -> VST, Wald padj) fires for an RNA-Seq
#   series with no expression values, or with `counts = TRUE` (GSE141295, GSE175759).
# #############################################################################
set.seed(1)
msg("Part C — per-dataset plots")

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) "p<0.001"
  else if (p < 0.01) paste0("p=", formatC(p, digits = 3, format = "f"))
  else paste0("p=", formatC(p, digits = 2, format = "f"))
}

ref_first <- function(grp, ref) c(ref, setdiff(names(sort(table(grp), decreasing = TRUE)), ref))

draw_gene_boxplot <- function(df, gene, co, out_dir, pfun, src) {
  lv  <- levels(df$Group)
  cmp <- lapply(setdiff(lv, co$ref), function(g) c(co$ref, g))
  cmp <- Filter(function(p) all(p %in% lv) &&
                  all(vapply(p, function(g) sum(df$Group == g) >= cfg$plots$min_n, TRUE)), cmp)
  ann <- vapply(cmp, function(p) fmt_p(pfun(df, p[1], p[2], gene)), character(1))
  rng <- range(df$Expression, na.rm = TRUE); span <- diff(rng)
  if (!is.finite(span) || span == 0) span <- 1
  ypos <- rng[2] + span * 0.10 * seq_along(cmp)

  p <- ggplot(df, aes(Group, Expression, fill = Group)) +
    geom_boxplot(outlier.shape = 16, outlier.size = 1.2, outlier.colour = "grey30",
                 outlier.alpha = 0.55, width = 0.6, linewidth = 0.55) +
    geom_jitter(width = 0.12, size = 2.3, alpha = 0.7) +
    scale_fill_manual(values = GROUP_PALETTE, na.value = "grey70") +
    labs(title = paste("Gene Expression –", gene),
         subtitle = sprintf("Dataset: %s (%s)", co$id, src),
         x = "Group", y = co$y_lab) +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
          legend.position = "none", plot.title = element_text(face = "bold"),
          panel.grid.major.x = element_blank())
  if (length(cmp))
    p <- p + ggsignif::geom_signif(comparisons = cmp, annotations = ann, y_position = ypos,
                                   tip_length = 0.01, textsize = 3.5, vjust = -0.2,
                                   colour = "grey35") +
             coord_cartesian(ylim = c(NA, max(ypos) + span * 0.10))
  ggsave(file.path(out_dir, paste0(co$id, "_boxplot_", gene, ".png")), p,
         width = cfg$plots$width, height = cfg$plots$height, dpi = cfg$dpi)
}

# Samples x genes heatmap, Z-score per gene, euclidean/complete clustering,
# samples annotated by group. Canvas sized to the content.
render_heatmap <- function(mat, grp, co, out_file, title, cellwidth = 11) {
  ord <- ref_first(grp, co$ref)
  ann <- data.frame(Group = factor(grp, levels = intersect(ord, unique(grp))),
                    row.names = rownames(mat))
  lv   <- levels(ann$Group)
  cols <- GROUP_PALETTE[lv]
  cols[is.na(cols)] <- grDevices::hcl.colors(sum(is.na(cols)), "Dark 3")
  names(cols) <- lv
  cellheight <- 5
  n_col <- 50                    # fixed colour steps, so 5 genes still get a gradient
  pheatmap::pheatmap(mat, scale = "column",
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",
                     clustering_method = "complete",
                     cellwidth = cellwidth, cellheight = cellheight,
                     color  = colorRampPalette(c("navy", "white", "red"))(n_col),
                     breaks = seq(-1.5, 1.5, length.out = n_col + 1),
                     annotation_row = ann, annotation_colors = list(Group = cols),
                     show_rownames = FALSE, show_colnames = TRUE,
                     fontsize = 14, fontsize_row = 11, fontsize_col = 11,
                     main = title, filename = out_file,
                     width  = if (ncol(mat) >= 20) 12 else ncol(mat) * cellwidth / 72 + 6,
                     height = nrow(mat) * cellheight / 72 + 7)
}

# Heatmap column labels: gene symbol where the platform maps one, else probe ID.
heatmap_labels <- function(ck, probes) {
  if (ck$co$idtype == "symbol") return(make.unique(probes))
  sym <- tryCatch(
    if (ck$co$idtype == "entrez")
      suppressMessages(AnnotationDbi::mapIds(getExportedValue("org.Hs.eg.db", "org.Hs.eg.db"),
                                             keys = sub("_at$", "", probes), column = "SYMBOL",
                                             keytype = "ENTREZID", multiVals = "first"))
    else ck$sym[probes],
    error = function(e) rep(NA_character_, length(probes)))
  make.unique(ifelse(is.na(sym), probes, unname(sym)))
}

draw_top_heatmap <- function(expr, grp, ck, out_dir, labeller) {
  expr <- expr[stats::complete.cases(expr), , drop = FALSE]     # hclust cannot take NA
  expr <- expr[rowSums(expr > median(expr)) >= 0.10 * ncol(expr), , drop = FALSE]
  n_top <- min(cfg$plots$heatmap_top, nrow(expr))
  top <- names(sort(apply(expr, 1, var), decreasing = TRUE))[seq_len(n_top)]
  mat <- t(expr[top, , drop = FALSE])
  colnames(mat) <- labeller(top)
  render_heatmap(mat, grp, ck$co, file.path(out_dir, paste0(ck$co$id, "_heatmap_top", n_top, ".png")),
                 sprintf("%s\nTop %d Variable Genes\n%s", ck$co$id, n_top, ck$co$y_lab))
}

draw_panel_heatmap <- function(panel, grp, co, out_dir) {
  panel <- panel[intersect(cfg$genes, rownames(panel)), , drop = FALSE]
  ok  <- stats::complete.cases(t(panel))
  render_heatmap(t(panel[, ok, drop = FALSE]), grp[ok], co,
                 file.path(out_dir, paste0(co$id, "_heatmap_", nrow(panel), "genes.png")),
                 sprintf("%s\nRenoRisk panel (%d genes)\n%s", co$id, nrow(panel), co$y_lab),
                 cellwidth = 40)
}

plot_gene_set <- function(mat, grp, ck, out_dir, pfun, src) {
  ord <- ref_first(grp, ck$co$ref)
  for (g in rownames(mat)) {
    df <- data.table(Expression = as.numeric(mat[g, ]),
                     Group = factor(grp, levels = intersect(ord, unique(grp))))
    draw_gene_boxplot(df[!is.na(Group) & !is.na(Expression)], g, ck$co, out_dir, pfun, src)
  }
}

plots_microarray <- function(ck, out_dir) {
  if (is.null(ck$panel)) stop("none of the panel genes map on this platform")
  labs_all <- group_labels(ck$sm$pheno, ck$co)
  grp <- labs_all[match(colnames(ck$panel), rownames(ck$sm$pheno))]
  wilcox_p <- function(d, a, b, ...) {
    x <- d$Expression[d$Group == a]; y <- d$Expression[d$Group == b]
    if (length(x) < cfg$plots$min_n || length(y) < cfg$plots$min_n) return(NA_real_)
    suppressWarnings(wilcox.test(y, x, exact = FALSE)$p.value)
  }
  plot_gene_set(ck$panel, grp, ck, out_dir, wilcox_p, unname(VENDOR[ck$co$idtype]))
  draw_panel_heatmap(ck$panel, grp, ck$co, out_dir)
  full <- on_log2(ck$sm$expr, ck$co)
  draw_top_heatmap(full, labs_all[match(colnames(full), rownames(ck$sm$pheno))], ck, out_dir,
                   function(pr) heatmap_labels(ck, pr))
  miss <- setdiff(cfg$genes, rownames(ck$panel))
  data.table(genes = nrow(ck$panel), samples = ncol(ck$panel), groups = uniqueN(grp),
             test = "Wilcoxon vs reference",
             missing_genes = if (length(miss)) paste(miss, collapse = ", ") else NA_character_)
}

# RNA-Seq: needs a counts file in the dataset folder (a series matrix carries
# no expression table for RNA-Seq). DESeq2 -> VST values, Wald padj vs reference.
plots_rnaseq <- function(ck, out_dir, folder) {
  cf <- list.files(folder, "count.*\\.(csv|tsv|txt)(\\.gz)?$", full.names = TRUE,
                   ignore.case = TRUE)
  if (!length(cf))
    stop("RNA-Seq cohort but no counts matrix in ", basename(folder),
         " - the series matrix has no expression table; supply counts or run a ",
         "quantification stage first")
  for (p in c("DESeq2", "SummarizedExperiment"))
    if (!requireNamespace(p, quietly = TRUE)) stop("package not installed: ", p)
  raw <- fread(cf[1])
  gid <- as.character(raw[[grep("gene|geneid|ensembl", names(raw), ignore.case = TRUE)[1]]])
  num <- names(raw)[vapply(raw, is.numeric, TRUE)]
  cnt <- round(as.matrix(raw[, ..num])); rownames(cnt) <- gid
  grp <- group_labels(ck$sm$pheno, ck$co)[match(colnames(cnt), rownames(ck$sm$pheno))]
  dds <- DESeq2::DESeqDataSetFromMatrix(cnt, data.frame(row.names = colnames(cnt),
                                                        condition = factor(grp)), ~ condition)
  dds <- DESeq2::DESeq(dds[rowSums(DESeq2::counts(dds)) >= 10, ], quiet = TRUE)
  vsd <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE))
  # VST shifts slightly with factor-level names/order; plot the analysis values instead.
  if (nrow(ck$sm$expr)) vsd <- ck$sm$expr[rownames(vsd), colnames(vsd), drop = FALSE]
  sym <- if (any(cfg$genes %in% rownames(vsd))) rownames(vsd) else
         suppressMessages(AnnotationDbi::mapIds(getExportedValue("org.Hs.eg.db", "org.Hs.eg.db"),
                                                rownames(vsd), "SYMBOL", "ENSEMBL",
                                                multiVals = "first"))
  hit <- which(sym %in% cfg$genes)
  if (!length(hit)) stop("no panel gene survived DESeq2 prefiltering")
  pnl <- vsd[hit, , drop = FALSE]; rownames(pnl) <- sym[hit]
  wald_p <- function(d, a, b, gene) {
    r  <- DESeq2::results(dds, contrast = c("condition", b, a))
    rn <- rownames(vsd)[hit][match(gene, sym[hit])]
    if (rn %in% rownames(r)) r[rn, "padj"] else NA_real_
  }
  plot_gene_set(pnl, grp, ck, out_dir, wald_p, "DNA-Seq")
  draw_top_heatmap(vsd, grp, ck, out_dir, identity)
  draw_panel_heatmap(pnl, grp, ck$co, out_dir)
  data.table(genes = length(hit), samples = ncol(vsd), groups = uniqueN(grp),
             test = "DESeq2 Wald padj")
}

plot_report <- list()
for (ck in CK) {
  out_dir <- file.path(cfg$out$plots, ck$co$dir)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  type <- if (ck$rnaseq) "rnaseq" else "microarray"
  # Normalised RNA-Seq values take the microarray path; DESeq2 only for counts.
  r <- tryCatch(if (ck$rnaseq && (!nrow(ck$sm$expr) || isTRUE(ck$co$counts)))
                  plots_rnaseq(ck, out_dir, file.path(cfg$raw_dir, ck$co$dir))
                else plots_microarray(ck, out_dir),
                error = function(e) { warning(ck$co$id, ": ", conditionMessage(e)); conditionMessage(e) })
  plot_report[[ck$co$id]] <- if (is.data.table(r))
    data.table(dataset = ck$co$dir, cohort = ck$co$id, type, status = "done", r)
  else data.table(dataset = ck$co$dir, cohort = ck$co$id, type, status = "failed",
                  genes = NA_integer_, samples = NA_integer_, groups = NA_integer_, test = r)
  msg(ck$co$dir, " [", type, "] -> ", out_dir)
}
write_xlsx(rbindlist(plot_report, fill = TRUE), file.path(cfg$out$plots, "boxplots_summary.xlsx"))
msg("per-dataset PNGs: ", length(list.files(cfg$out$plots, "\\.png$", recursive = TRUE)))
write_session(cfg$out$plots, "boxplots")


# #############################################################################
# PART D — STAGE 3: BATCH CORRECTION (and the Stage 2 hand-off tables)
#   One limma::removeBatchEffect model over the aggregated data: every dataset
#   except the subsets of GSE30122 (SUBSETS), so no sample enters twice.
#     batch  = StudyID (17 studies; plan §3 "~ StudyID + Condition")
#     design = intercept + Disease
#   The per-dataset tables (cohort_qc, sample_indication) and the uncorrected
#   hand-off (panel_uncorrected_long, input of the per-study effects) keep all
#   14 datasets; the corrected values and every figure here are aggregated.
#   Compartment main effects are aliased with batch and NOT in the design;
#   removeBatchEffect would silently zero aliased batch coefficients, so a rank
#   check guards it. Fitted gene by gene on the samples that have the gene
#   (C4BPA is absent from GSE142025). Corrects location, not spread.
# #############################################################################
set.seed(1)
msg("Part D — Stage 2/3")

panels <- list(); phenos <- list(); qc <- list(); indic <- list()
for (ck in CK) {
  if (is.null(ck$panel)) next
  co <- ck$co; mat <- ck$panel
  ph <- data.table(GSM = colnames(mat), StudyID = co$id,
                   Condition   = derive_condition(ck$sm$pheno, co),
                   Compartment = derive_compartment(ck$sm$pheno, co))
  ph <- cbind(ph, derive_indication(ck$sm$pheno, co, ph$Condition))[!is.na(Compartment)]
  mat <- mat[, ph$GSM, drop = FALSE]
  indic[[co$id]]  <- ph[, .(StudyID, GSM, Compartment, Condition, Indication, Label)]
  phenos[[co$id]] <- ph[, .(GSM, StudyID, Condition, Compartment)]
  panels[[co$id]] <- mat
  qc[[co$id]] <- data.table(StudyID = co$id, genes = nrow(mat), n = ncol(mat),
                            Control = sum(ph$Condition == "Control"),
                            Disease = sum(ph$Condition == "Disease"),
                            median  = round(median(mat, na.rm = TRUE), 2))
}
stopifnot(length(panels) > 1)
write_xlsx(rbindlist(qc), file.path(cfg$out$s23, "cohort_qc.xlsx"))
indic <- rbindlist(indic)
write_xlsx(indic, file.path(cfg$out$s23, "sample_indication.xlsx"))
msg("indications: ", paste(indic[, .N, by = Indication][order(-N), sprintf("%s %d", Indication, N)],
                           collapse = ", "))

# A gene missing from one cohort (C4BPA in the GSE142025 RNA-Seq deposit) stays
# in the panel: that cohort's values are NA and every per-gene step below uses
# the samples that have the gene.
genes  <- intersect(unique(unlist(lapply(panels, rownames))), cfg$genes)   # row order as before
common <- Reduce(intersect, lapply(panels, rownames))   # complete genes: PCA / UMAP
msg("panel genes: ", paste(genes, collapse = ", "), "; on every cohort: ",
    paste(intersect(cfg$genes, common), collapse = ", "))

pheno <- rbindlist(phenos)
expr  <- do.call(cbind, lapply(panels, function(m) {
  out <- matrix(NA_real_, length(genes), ncol(m), dimnames = list(genes, colnames(m)))
  hit <- intersect(genes, rownames(m)); out[hit, ] <- m[hit, , drop = FALSE]; out
}))
# GSE30122 shares its GSM accessions with GSE30528/GSE30529: column names are
# not unique, so samples are indexed by position, never by name.
stopifnot(identical(unname(colnames(expr)), pheno$GSM))
pheno[, col := .I]
pheno[, Batch := if (uniqueN(Compartment) > 1) paste(StudyID, Compartment, sep = ":")
                 else StudyID, by = StudyID]
ok <- pheno[, .N, by = .(Batch, Condition)][N >= cfg$s23$min_group_n]
ph <- pheno[paste(Batch, Condition) %in% ok[, paste(Batch, Condition)]]
if (nrow(ph) < nrow(pheno))
  msg("dropped ", nrow(pheno) - nrow(ph), " samples in arms with n < ", cfg$s23$min_group_n)
stopifnot(uniqueN(ph$Batch) > 1, uniqueN(ph$Condition) == 2)

# Hand-off to Parts E–H: every dataset's samples BEFORE correction (effects are
# modelled or computed within study), all identifier columns as text. Each
# later part drops SUBSETS itself wherever it aggregates.
to_long <- function(d, m, ids)
  melt(cbind(d[, ..ids], t(m)), id.vars = ids, variable.name = "Gene",
       value.name = "Expression")[!is.na(Expression)]
unc <- to_long(ph, expr[, ph$col, drop = FALSE], c("GSM", "StudyID", "Condition", "Compartment", "Batch"))
unc[, `:=`(Condition = as.character(Condition), Gene = as.character(Gene))]
write_xlsx(unc, file.path(cfg$out$s23, "panel_uncorrected_long.xlsx"))

# The aggregated data: SUBSETS left out (their samples are already in GSE30122).
pph <- ph[!StudyID %in% SUBSETS]
pph[, Compartment := factor(Compartment)]
batch <- factor(pph$StudyID)   # plan §3 Stage 3: ~ StudyID + Condition
# The same model, fitted gene by gene on the samples that have the gene: a batch
# unit with no values for a gene would otherwise make the design rank-deficient
# and removeBatchEffect would silently zero a batch coefficient. For a gene
# present everywhere this equals the all-genes call (it fits genes independently).
adj <- expr[, pph$col, drop = FALSE]
for (g in rownames(adj)) {
  ok <- !is.na(adj[g, ])
  b  <- droplevels(batch[ok])
  des <- model.matrix(~ Condition, data = pph[ok])
  X.batch <- model.matrix(~ b, contrasts.arg = list(b = contr.sum(levels(b))))[, -1, drop = FALSE]
  stopifnot("design + batch is rank-deficient: some batch effects would be silently zeroed" =
              qr(cbind(des, X.batch))$rank == ncol(des) + ncol(X.batch))
  adj[g, ok] <- limma::removeBatchEffect(adj[g, ok, drop = FALSE], batch = b, design = des)
}
msg("one model: corrected ", ncol(adj), " samples across ", nlevels(batch),
    " batch units = studies (", uniqueN(pph$StudyID), " studies, ", nlevels(pph$Compartment),
    " compartments)", if (length(SUBSETS)) paste0("; not aggregated: ", paste(SUBSETS, collapse = ", ")))
pph[, Compartment := as.character(Compartment)]
stopifnot(!anyDuplicated(pph$GSM))       # no sample counted twice

harm <- to_long(pph, adj, c("GSM", "StudyID", "Condition", "Compartment"))
harm[, Gene := factor(Gene, levels = genes)]
write_xlsx(harm, file.path(cfg$out$s23, "panel_harmonised_long.xlsx"))

# Residual spread per study, per gene (pooling genes would look falsely uniform).
sd_tab <- harm[, .(sd = sd(Expression, na.rm = TRUE)), by = .(Compartment, Gene, StudyID)]
sd_rng <- sd_tab[, .(min_sd = round(min(sd), 2), max_sd = round(max(sd), 2),
                     ratio = round(max(sd) / min(sd), 1)), by = .(Compartment, Gene)][order(-ratio)]
write_xlsx(sd_tab, file.path(cfg$out$s23, "residual_sd_by_study_gene.xlsx"))
write_xlsx(sd_rng, file.path(cfg$out$s23, "residual_sd_ratio.xlsx"))
msg("residual SD ratio across studies after correction: max ", max(sd_rng$ratio))

# --- harmonised panel figures -------------------------------------------------
sub_lab <- sprintf("%d cohorts, %d samples", uniqueN(harm$StudyID), uniqueN(harm$GSM))
n_units <- nlevels(batch)
grid_theme <- function() {
  if (cfg$s23$gridlines %in% c("horizontal", "both"))
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
  else theme(panel.grid = element_blank())
}
# Separators between groups are vlines on the boundaries (x = 1.5, 2.5, ...);
# discrete-axis gridlines would cut through the middle of each pair.
sep_layer <- function(d)
  if (cfg$s23$gridlines %in% c("vertical", "both"))
    geom_vline(data = d, aes(xintercept = xintercept), inherit.aes = FALSE,
               colour = "grey78", linewidth = 0.3)
sep_gene  <- data.table(xintercept = seq_len(nlevels(harm$Gene) - 1L) + 0.5)
sep_study <- harm[, .(n = uniqueN(StudyID)), by = Compartment][
  n > 1, .(xintercept = seq_len(n - 1L) + 0.5), by = Compartment]

p <- ggplot(harm, aes(Gene, Expression, fill = Condition)) +
  sep_layer(sep_gene) +
  geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.45,
               position = position_dodge(0.78), alpha = 0.9) +
  # thin dense strata so the points do not swamp the boxes
  geom_point(data = harm[, .SD[if (.N > 120) sample(.N, 120) else seq_len(.N)],
                         by = .(Compartment, Gene, Condition)],
             position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.78),
             size = 0.55, alpha = 0.22, colour = "grey20") +
  facet_wrap(~ Compartment, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = COND_PALETTE) +
  labs(title = "RenoRisk panel across curated GEO kidney-biopsy cohorts", subtitle = sub_lab,
       x = NULL, y = "Expression (log2)",
       caption = "Submitter-processed series-matrix values; study offsets removed, scale differences not.") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 9),
        plot.caption = element_text(colour = "grey45", size = 7.5),
        axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
        axis.title.y = element_text(size = 10), legend.position = "top") +
  grid_theme()
ggsave(file.path(cfg$out$s23, "panel_boxplot_by_compartment.png"), p,
       width = cfg$s23$width, height = cfg$s23$height, dpi = cfg$dpi)

# Per-cohort view: did the correction align the studies?
p <- ggplot(harm, aes(StudyID, Expression, fill = Condition)) +
  sep_layer(sep_study) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4, linewidth = 0.4,
               width = 0.36, position = position_dodge(0.75)) +
  facet_grid(Gene ~ Compartment, scales = "free", space = "free_x") +
  scale_fill_manual(values = COND_PALETTE) +
  labs(title = "Per-cohort check after batch correction", subtitle = sub_lab,
       x = NULL, y = "Batch-corrected expression (log2)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(5.5, 5.5, 5.5, 34), legend.position = "top") +
  grid_theme()
ggsave(file.path(cfg$out$s23, "panel_boxplot_per_cohort.png"), p,
       width = cfg$s23$pc_width, height = cfg$s23$pc_height, dpi = cfg$dpi)

# --- Plot 1: PCA and UMAP before vs after correction ---------------------------
# Same samples and genes, each embedding fitted separately per matrix; genes
# scaled to unit variance (per-gene-centred cohorts make raw variances
# incomparable). Batch removal is quantified: R² of batch / Condition on PC1-2,
# and the mean silhouette width of the batch units on the UMAP.
pca_one <- function(m, when) {
  pc <- prcomp(t(m), center = TRUE, scale. = TRUE)
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  sc <- data.table(pph[, .(StudyID, Batch, Compartment, Condition)],
                   PC1 = pc$x[, 1], PC2 = pc$x[, 2], When = when)
  r2 <- rbindlist(lapply(1:2, function(k) data.table(
    When = when, PC = paste0("PC", k), var_explained_pct = round(ve[k], 1),
    R2_batch     = round(summary(lm(sc[[paste0("PC", k)]] ~ sc$Batch))$r.squared, 3),
    R2_condition = round(summary(lm(sc[[paste0("PC", k)]] ~ sc$Condition))$r.squared, 3))))
  list(scores = sc, r2 = r2)
}
umap_one <- function(m, when) {
  xy  <- uwot::umap(scale(t(m)), n_neighbors = cfg$s23$umap_neighbors,
                    min_dist = cfg$s23$umap_min_dist, n_threads = 1, seed = 1)
  sil <- mean(cluster::silhouette(as.integer(factor(pph$Batch)), dist(xy))[, "sil_width"])
  list(scores = data.table(pph[, .(StudyID, Batch, Compartment, Condition)],
                           UMAP1 = xy[, 1], UMAP2 = xy[, 2], When = when),
       sil = data.table(When = when, batch_silhouette = round(sil, 3)))
}
embed_plot <- function(scores, x, y, title, caption, file) {
  p <- ggplot(scores, aes(.data[[x]], .data[[y]], colour = StudyID, shape = Compartment)) +
    geom_point(size = 1.3, alpha = 0.75) +
    facet_wrap(~ When, nrow = 1, scales = "free") +
    scale_shape_manual(values = c(16, 17, 15)) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1), ncol = 1),
           shape  = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    labs(title = title, subtitle = sub_lab, caption = caption, x = x, y = y) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey35", size = 9),
          plot.caption = element_text(colour = "grey45", size = 7.5),
          strip.text = element_text(size = 8.5))
  ggsave(file.path(cfg$out$s23, file), p, width = cfg$s23$pca_width,
         height = cfg$s23$pca_height, dpi = cfg$dpi)
}
relabel <- function(scores, labs) {
  scores[labs, When := i.lab, on = "When"]
  scores[, When := factor(When, levels = labs$lab)]
}

cg  <- intersect(cfg$genes, common)   # PCA needs complete rows: genes on every cohort
pcs <- list(pca_one(expr[cg, pph$col, drop = FALSE], "Before correction"),
            pca_one(adj[cg, , drop = FALSE], "After correction"))
pca_r2 <- rbindlist(lapply(pcs, `[[`, "r2"))
write_xlsx(pca_r2, file.path(cfg$out$s23, "pca_batch_r2.xlsx"))
scores <- relabel(rbindlist(lapply(pcs, `[[`, "scores")),
  pca_r2[, .(lab = sprintf("%s\nPC1 %.0f%%, PC2 %.0f%% of variance; batch R² PC1 %.2f, PC2 %.2f",
                           When[1], var_explained_pct[1], var_explained_pct[2],
                           R2_batch[1], R2_batch[2])), by = When])
embed_plot(scores, "PC1", "PC2",
           sprintf("PCA of the panel (%s) before and after batch correction", paste(cg, collapse = ", ")),
           paste("Genes centred and scaled to unit variance; PCA fitted separately per panel.",
                 sprintf("Batch R² = share of the PC explained by the %d study × compartment units.",
                         n_units)),
           "pca_before_after_correction.png")
msg("PCA batch R² on PC1: ", pca_r2[PC == "PC1", paste(When, R2_batch, collapse = " -> ")])

ums <- list(umap_one(expr[cg, pph$col, drop = FALSE], "Before correction"),
            umap_one(adj[cg, , drop = FALSE], "After correction"))
umap_sil <- rbindlist(lapply(ums, `[[`, "sil"))
write_xlsx(umap_sil, file.path(cfg$out$s23, "umap_batch_silhouette.xlsx"))
scores <- relabel(rbindlist(lapply(ums, `[[`, "scores")),
  umap_sil[, .(lab = sprintf("%s\nbatch silhouette %.2f", When, batch_silhouette)), by = When])
embed_plot(scores, "UMAP1", "UMAP2",
           sprintf("UMAP of the panel (%s) before and after batch correction", paste(cg, collapse = ", ")),
           sprintf(paste("Genes centred and scaled to unit variance; UMAP fitted separately per panel",
                         "(n_neighbors = %d, min_dist = %.2f, seed 1). Batch silhouette = mean",
                         "silhouette width of the %d study × compartment units."),
                   cfg$s23$umap_neighbors, cfg$s23$umap_min_dist, n_units),
           "umap_before_after_correction.png")
write_session(cfg$out$s23, "stage23")


# #############################################################################
# PART E — STAGE 4: MULTI-FACTOR STATISTICAL MODELLING
#   Per gene: Expression ~ 0 + Batch + Disease:Compartment + Sex, two-step WLS
#   (OLS -> residual variance per batch unit -> weights 1/variance), on the
#   UNCORRECTED values (batch fitted, not pre-removed). Reports the Disease
#   effect per compartment, the between-compartment differences of that effect
#   and sex; BH within each family.
#   Sex: reported (GSE116626, GSE129973, GSE96804) or inferred from XIST vs Y
#   genes; inference accepted only for clearly bimodal cohorts.
#   Single-dataset calculations (sex per cohort, per-unit effects) cover all 14
#   datasets (GSE142025 RNA-Seq included); the multi-factor model aggregates,
#   so it leaves SUBSETS out (1078 samples, 18 units).
# #############################################################################
set.seed(1)
msg("Part E — Stage 4")
samples <- unique(unc[, .(GSM, StudyID, Condition, Compartment, Batch)])

reported_sex <- function(pheno) {
  v <- tryCatch(field_values(pheno, "sex"), error = function(e) NULL)
  if (is.null(v)) return(rep(NA_character_, nrow(pheno)))
  v <- tolower(trimws(v))
  ifelse(v %in% c("m", "male"), "M", ifelse(v %in% c("f", "female"), "F", NA_character_))
}

# Score = mean z(Y genes) - z(XIST) within cohort (whichever side the platform
# carries), split by 2-means; accepted only if the clusters are well separated.
infer_sex <- function(ck, gsm) {
  map <- gene_probes(ck, SEX_GENES)
  y <- intersect(cfg$s4$y_genes, map); has_x <- cfg$s4$x_gene %in% map
  if (!has_x && !length(y))
    return(list(sex = rep(NA_character_, length(gsm)), sep = NA_real_,
                note = "sex marker genes not on platform"))
  ex <- on_log2(ck$sm$expr, ck$co)
  z  <- t(scale(t(collapse_genes(ex[, colnames(ex) %in% gsm, drop = FALSE], map))))
  score <- (if (length(y)) colMeans(z[y, , drop = FALSE]) else 0) -
           (if (has_x) z[cfg$s4$x_gene, ] else 0)
  km  <- kmeans(score, centers = range(score))
  sep <- abs(diff(km$centers[, 1])) /
         sqrt(mean(tapply(score, km$cluster, function(s) if (length(s) > 1) var(s) else 0)))
  if (!is.finite(sep) || sep < cfg$s4$min_sex_sep)
    return(list(sex = setNames(rep(NA_character_, length(score)), names(score)), sep = sep,
                note = "no clear XIST/Y bimodality"))
  list(sex = setNames(ifelse(km$cluster == which.max(km$centers[, 1]), "M", "F"), names(score)),
       sep = sep, note = paste(c(if (length(y)) sprintf("%d Y genes", length(y)),
                                 if (has_x) "XIST"), collapse = " + "))
}

sex_tab <- list(); sex_log <- list()
for (sid in unique(samples$StudyID)) {
  ck  <- CK[[sid]]
  gsm <- colnames(ck$sm$expr)[colnames(ck$sm$expr) %in% samples[StudyID == sid, GSM]]
  rep_sex <- reported_sex(ck$sm$pheno)
  if (isTRUE(ck$co$title_sex))   # e.g. GSE93798 "... H12_M_67_1_Con_L"
    rep_sex <- sub("^.* [A-Z]+[0-9]+_([MF])_.*$", "\\1", ck$sm$pheno$title)
  rep_sex <- setNames(rep_sex, ck$sm$pheno$geo_accession)[gsm]
  stopifnot(all(rep_sex %in% c("M", "F", NA)))
  inf <- infer_sex(ck, gsm)
  sex_tab[[sid]] <- data.table(StudyID = sid, GSM = gsm, Reported = unname(rep_sex),
                               Inferred = unname(inf$sex[gsm]))
  sex_log[[sid]] <- data.table(StudyID = sid, n = length(gsm), separation = round(inf$sep, 1),
                               inference = inf$note, reported = sum(!is.na(rep_sex)))
}
sex_tab <- rbindlist(sex_tab)
sex_tab[, Sex := fifelse(!is.na(Reported), Reported, Inferred)]
sex_tab[, SexSource := fifelse(!is.na(Reported), "reported",
                       fifelse(!is.na(Inferred), "inferred", "unknown"))]
conc <- sex_tab[!is.na(Reported) & !is.na(Inferred),
                .(agree = sum(Reported == Inferred),
                  concordance = round(mean(Reported == Inferred), 3)), by = StudyID]
sex_log <- merge(rbindlist(sex_log), conc, by = "StudyID", all.x = TRUE)
sex_log <- merge(sex_log, sex_tab[, .(male = sum(Sex == "M", na.rm = TRUE),
                                      female = sum(Sex == "F", na.rm = TRUE),
                                      unknown = sum(is.na(Sex))), by = StudyID], by = "StudyID")
write_xlsx(sex_log, file.path(cfg$out$s4, "sex_by_cohort.xlsx"))
write_xlsx(sex_tab, file.path(cfg$out$s4, "sex_by_sample.xlsx"))

dat <- merge(unc, sex_tab[, .(StudyID, GSM, Sex, SexSource)], by = c("StudyID", "GSM"))
n_unknown <- uniqueN(dat[is.na(Sex), .(StudyID, GSM)])
if (n_unknown) msg("excluded, sex unknown: ", n_unknown, " observations")
dat <- dat[!is.na(Sex)]
dat[, Gene := factor(Gene, levels = intersect(cfg$genes, unique(Gene)))]
mdat <- dat[!StudyID %in% SUBSETS]       # the aggregated model's data
stopifnot(!anyDuplicated(mdat[Gene == levels(Gene)[1], GSM]))
msg("model: ", uniqueN(mdat$GSM), " samples, ", uniqueN(mdat$Batch), " batch units, ",
    nlevels(mdat$Gene), " genes")
comps <- sort(unique(mdat$Compartment))

build_design <- function(d) {
  X <- model.matrix(~ 0 + factor(Batch), data = d)
  colnames(X) <- sub("^factor\\(Batch\\)", "Batch:", colnames(X))
  for (cp in comps)
    X <- cbind(X, matrix(as.numeric(d$Condition == "Disease" & d$Compartment == cp),
                         dimnames = list(NULL, paste0("Disease:", cp))))
  if (uniqueN(d$Sex) > 1) X <- cbind(X, `Sex:M` = as.numeric(d$Sex == "M"))
  X
}

# Disease effect per compartment, sex, and every pairwise difference of the
# Disease effects, as linear contrasts of the coefficients.
contrasts_of <- function(est, V, df) {
  terms <- grep("^(Disease|Sex):", names(est), value = TRUE)
  L <- diag(length(est))[match(terms, names(est)), , drop = FALSE]
  rownames(L) <- terms
  dterms <- grep("^Disease:", names(est), value = TRUE)
  for (a in seq_along(dterms)) for (b in seq_along(dterms)) if (a < b) {
    l <- setNames(numeric(length(est)), names(est))
    l[dterms[b]] <- 1; l[dterms[a]] <- -1
    L <- rbind(L, l)
    rownames(L)[nrow(L)] <- sprintf("Diff:%s-%s", sub("Disease:", "", dterms[b]),
                                    sub("Disease:", "", dterms[a]))
  }
  e <- drop(L %*% est); se <- sqrt(diag(L %*% V %*% t(L))); t <- e / se
  data.table(Term = rownames(L), log2FC = e, SE = se, t = t, df = df, p = 2 * pt(-abs(t), df))
}

fit_gene <- function(d) {
  X <- build_design(d)
  stopifnot("design is rank-deficient" = qr(X)$rank == ncol(X))
  y  <- d$Expression
  r0 <- lm.fit(X, y)$residuals
  v  <- tapply(r0^2, d$Batch, function(r) sum(r) / max(length(r) - 1, 1))
  fit <- lm(y ~ 0 + X, weights = 1 / v[d$Batch])
  est <- setNames(coef(fit), colnames(X))
  V   <- vcov(fit); dimnames(V) <- list(colnames(X), colnames(X))
  contrasts_of(est, V, fit$df.residual)
}

res <- rbindlist(lapply(levels(mdat$Gene), function(g) cbind(Gene = g, fit_gene(mdat[Gene == g]))))
res[, Family := sub(":.*", "", Term)]
res[, FDR := p.adjust(p, "BH"), by = Family]
res[, Term := sub("^Disease:", "Disease vs Control | ", Term)]
res[, Term := sub("^Diff:", "Disease effect difference | ", Term)]
res[, Term := sub("^Sex:M", "Sex (M vs F)", Term)]
res <- res[, .(Gene, Family, Term, log2FC = signif(log2FC, 4), SE = signif(SE, 4),
               t = round(t, 3), df, p = signif(p, 4), FDR = signif(FDR, 4))]
write_xlsx(res, file.path(cfg$out$s4, "stage4_multifactor_results.xlsx"))
msg("Stage 4 disease effects with FDR < 0.05: ", res[Family == "Disease" & FDR < 0.05, .N], " of ",
    res[Family == "Disease", .N])

# Per-unit Disease effects (sex-adjusted where sex varies): one dataset at a
# time, so all 14 datasets.
per_study <- dat[, {
  if (uniqueN(Condition) != 2 || min(table(Condition)) < 2) NULL else {
    adj_sex <- uniqueN(Sex) > 1 && all(table(Sex) >= 2)
    s <- summary(if (adj_sex) lm(Expression ~ Condition + Sex)
                 else lm(Expression ~ Condition))$coefficients["ConditionDisease", ]
    .(n_control = sum(Condition == "Control"), n_disease = sum(Condition == "Disease"),
      sex_adjusted = adj_sex, log2FC = signif(s[1], 4), SE = signif(s[2], 4),
      p = signif(s[4], 4))
  }
}, by = .(Gene, Batch, Compartment)]
write_xlsx(per_study, file.path(cfg$out$s4, "stage4_per_study_effects.xlsx"))

pd <- res[Family == "Disease"]
pd[, `:=`(Compartment = sub("Disease vs Control \\| ", "", Term),
          lo = log2FC - qt(0.975, df) * SE, hi = log2FC + qt(0.975, df) * SE,
          Gene = factor(Gene, levels = rev(cfg$genes)),
          sig = fifelse(FDR < 0.05, "FDR < 0.05", "n.s."))]
p <- ggplot(pd, aes(log2FC, Gene, colour = sig)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, orientation = "y") +
  geom_point(size = 2.4) +
  facet_wrap(~ Compartment, nrow = 1) +
  scale_colour_manual(values = c(`FDR < 0.05` = "#B2182B", n.s. = "grey45"), name = NULL) +
  labs(title = "Stage 4 — Disease vs Control, multi-factor model",
       subtitle = sprintf("Expression ~ batch unit + Disease per compartment + sex (WLS); %d samples, %d batch units",
                          uniqueN(mdat$GSM), uniqueN(mdat$Batch)),
       x = "log2 fold change (95% CI)", y = NULL,
       caption = "Submitter-processed values; sex reported for 3 cohorts, inferred from XIST / Y genes otherwise.") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top",
        plot.subtitle = element_text(colour = "grey35", size = 9),
        plot.caption = element_text(colour = "grey45", size = 7.5))
ggsave(file.path(cfg$out$s4, "stage4_disease_effects.png"), p,
       width = cfg$s4$width, height = cfg$s4$height, dpi = cfg$dpi)
write_session(cfg$out$s4, "stage4")


# #############################################################################
# PART F — POOLING STRATEGY (plan §5): RANDOM-EFFECTS META-ANALYSIS
#   Unit = batch unit. Per unit and gene: Hedges' g (bias-corrected SMD,
#   primary) and log2FC (MD, secondary) with variances (metafor::escalc) —
#   computed for all 15 units (14 microarray + GSE142025 RNA-Seq; single-dataset
#   calculation; column `Pooled`). Pooled by DerSimonian–Laird (metafor::rma)
#   over all 18 units (Tubulointerstitial k = 4, Whole biopsy k = 4, 3 for C4BPA):
#   all units, per compartment, compartment moderator; Q, I², tau², 95%
#   prediction interval; leave-one-out. Same rule for the indication and
#   control-tissue pooling.
#   Figures: forest plot per gene (Plot 2), summary forest, indication heatmap
#   (Plot 4), control-tissue sensitivity. Input = `unc` (within-study effects
#   need no batch correction).
# #############################################################################
set.seed(1)
msg("Part F — pooling strategy")
pmin_n <- cfg$pool$min_n

# Disease vs Control effect sizes from per-arm summaries (`m_`, `s_`, `n_`).
effect_sizes <- function(w, measure)
  escalc(measure = measure, m1i = m_Disease, sd1i = s_Disease, n1i = n_Disease,
         m2i = m_Control, sd2i = s_Control, n2i = n_Control, data = w)
arm_summaries <- function(d, by)
  dcast(d[, .(m = mean(Expression), s = sd(Expression), n = .N), by = c(by, "Condition")],
        as.formula(paste(paste(by, collapse = " + "), "~ Condition")),
        value.var = c("m", "s", "n"))[n_Control >= pmin_n & n_Disease >= pmin_n &
                                        s_Control > 0 & s_Disease > 0]
pool_dl <- function(yi, vi, ...) rma(yi, vi, method = cfg$pool$method, ...)
wald_ci <- function(y, v) list(lb = y - qnorm(0.975) * sqrt(v), ub = y + qnorm(0.975) * sqrt(v))
sig4 <- function(d) d[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 4) else x)]

# --- 1. per-unit effect sizes ---------------------------------------------------
wide <- arm_summaries(unc, c("Gene", "Batch", "StudyID", "Compartment"))
es_g <- effect_sizes(wide, "SMD"); es_fc <- effect_sizes(wide, "MD")
eff <- data.table(wide[, .(Gene, Batch, StudyID, Compartment, n_Control, n_Disease)],
                  g = as.numeric(es_g$yi), g_var = as.numeric(es_g$vi),
                  log2FC = as.numeric(es_fc$yi), log2FC_var = as.numeric(es_fc$vi))
eff[, `:=`(Gene = factor(Gene, levels = intersect(cfg$genes, unique(Gene))),
           Pooled = !StudyID %in% SUBSETS)]
setorder(eff, Gene, Compartment, Batch)
write_xlsx(eff, file.path(cfg$out$pool, "pooling_study_effects.xlsx"))
msg(nrow(eff), " gene x unit effect sizes from ", uniqueN(eff$Batch), " units; pooled: ",
    uniqueN(eff[(Pooled), Batch]), " units")

# --- 2. DL pooling: all units, per compartment, moderator ----------------------
summarise_rma <- function(fit, gene, measure, subset) {
  pr <- predict(fit)
  data.table(Gene = gene, Measure = measure, Subset = subset, k = fit$k,
             estimate = as.numeric(fit$beta), SE = fit$se, ci_lb = fit$ci.lb, ci_ub = fit$ci.ub,
             pi_lb = pr$pi.lb, pi_ub = pr$pi.ub, z = fit$zval, p = fit$pval,
             Q = fit$QE, Q_df = fit$k - 1, Q_p = fit$QEp, I2 = fit$I2, tau2 = fit$tau2)
}
pooled <- list(); moder <- list(); fits <- list()
for (gene in levels(eff$Gene)) for (ms in c("g", "log2FC")) {
  d  <- eff[Gene == gene & Pooled]
  yi <- d[[ms]]; vi <- d[[paste0(ms, "_var")]]
  fit <- pool_dl(yi, vi)
  if (ms == "g") fits[[gene]] <- list(fit = fit, d = d)
  pooled[[length(pooled) + 1]] <- summarise_rma(fit, gene, ms, "All units")
  for (cp in sort(unique(d$Compartment))) {
    i <- d$Compartment == cp
    pooled[[length(pooled) + 1]] <- if (sum(i) >= 2)
      summarise_rma(pool_dl(yi[i], vi[i]), gene, ms, cp)
    else with(wald_ci(yi[i], vi[i]),          # one unit: the study estimate itself
      data.table(Gene = gene, Measure = ms, Subset = cp, k = 1L, estimate = yi[i],
                 SE = sqrt(vi[i]), ci_lb = lb, ci_ub = ub, z = yi[i] / sqrt(vi[i]),
                 p = 2 * pnorm(-abs(yi[i] / sqrt(vi[i])))))
  }
  if (uniqueN(d$Compartment) > 1) {
    mf <- pool_dl(yi, vi, mods = ~ Compartment, data = d)
    moder[[length(moder) + 1]] <- data.table(Gene = gene, Measure = ms, k = mf$k, QM = mf$QM,
      QM_df = mf$QMdf[1], QM_p = mf$QMp, residual_I2 = mf$I2)
  }
}
pooled <- rbindlist(pooled, fill = TRUE)
pooled[, FDR := p.adjust(p, "BH"), by = .(Measure, Subset)]     # across the 5 genes
moder <- rbindlist(moder)
moder[, FDR := p.adjust(QM_p, "BH"), by = Measure]
num <- names(pooled)[vapply(pooled, is.numeric, TRUE) & !names(pooled) %in% c("k", "Q_df")]
write_xlsx(copy(pooled)[, (num) := lapply(.SD, signif, 4), .SDcols = num],
           file.path(cfg$out$pool, "pooling_estimates.xlsx"))
write_xlsx(sig4(moder), file.path(cfg$out$pool, "pooling_compartment_moderator.xlsx"))
msg("pooled Hedges' g, all units: ",
    pooled[Measure == "g" & Subset == "All units",
           paste(sprintf("%s %.2f (FDR %.3g, I² %.0f%%)", Gene, estimate, FDR, I2), collapse = "; ")])

# --- 3. leave-one-out (all units, Hedges' g) ------------------------------------
loo <- rbindlist(lapply(names(fits), function(gene) {
  l <- leave1out(fits[[gene]]$fit)
  data.table(Gene = gene, Omitted = fits[[gene]]$d$Batch, g = l$estimate,
             ci_lb = l$ci.lb, ci_ub = l$ci.ub, p = l$pval, I2 = l$I2, tau2 = l$tau2)
}))
write_xlsx(sig4(loo), file.path(cfg$out$pool, "pooling_leave1out.xlsx"))

# --- 4. forest plot per gene (Plot 2): units by compartment, subgroup and
#        overall diamonds, overall on top ---------------------------------------
diamond <- function(est, lb, ub, y, h = 0.35)
  data.table(x = c(lb, est, ub, est), y = c(y, y + h, y, y - h))
ci_txt <- function(e, l, u) sprintf("%.2f [%.2f, %.2f]", e, l, u)

for (gene in levels(eff$Gene)) {
  d  <- eff[Gene == gene & Pooled]            # the plot shows what is pooled
  d[, c("lo", "hi") := wald_ci(g, g_var)]
  pg <- pooled[Gene == gene & Measure == "g"]
  rows <- list(); dias <- list(); y <- 0
  for (cp in rev(sort(unique(d$Compartment)))) {
    sg <- pg[Subset == cp]; y <- y + 1
    dias[[length(dias) + 1]] <- cbind(diamond(sg$estimate, sg$ci_lb, sg$ci_ub, y),
                                      grp = cp, kind = "Subgroup")
    rows[[length(rows) + 1]] <- data.table(y = y, label = sprintf("%s (k=%d)", cp, sg$k),
      est = sg$estimate, lo = sg$ci_lb, hi = sg$ci_ub, kind = "Subgroup", n = "",
      txt = ci_txt(sg$estimate, sg$ci_lb, sg$ci_ub))
    dd <- d[Compartment == cp][order(-g)]
    for (i in seq_len(nrow(dd))) {
      y <- y + 1
      rows[[length(rows) + 1]] <- data.table(y = y, label = paste0("   ", dd$Batch[i]),
        est = dd$g[i], lo = dd$lo[i], hi = dd$hi[i], kind = "Study",
        n = sprintf("%d / %d", dd$n_Control[i], dd$n_Disease[i]),
        txt = ci_txt(dd$g[i], dd$lo[i], dd$hi[i]))
    }
    y <- y + 0.6
  }
  ov <- pg[Subset == "All units"]; y <- y + 1
  dias[[length(dias) + 1]] <- cbind(diamond(ov$estimate, ov$ci_lb, ov$ci_ub, y, 0.45),
                                    grp = "All", kind = "Overall")
  rows[[length(rows) + 1]] <- data.table(y = y, label = sprintf("Overall RE model (k=%d)", ov$k),
    est = ov$estimate, lo = ov$ci_lb, hi = ov$ci_ub, kind = "Overall", n = "",
    txt = ci_txt(ov$estimate, ov$ci_lb, ov$ci_ub))
  rows <- rbindlist(rows); dias <- rbindlist(dias)
  ymax <- max(rows$y)
  rows[, y := ymax + 1 - y]; dias[, y := ymax + 1 - y]
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
    coord_cartesian(xlim = range(c(rows$lo, rows$hi, 0), na.rm = TRUE)) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          axis.ticks.y = element_blank(), plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey35", size = 9),
          plot.caption = element_text(colour = "grey45", size = 7),
          axis.title.y.right = element_text(size = 9, angle = 0, vjust = 1))
  ggsave(file.path(cfg$out$pool, sprintf("pooling_forest_%s.png", gene)), p,
         width = cfg$pool$width, height = 1.6 + 0.27 * nrow(rows), dpi = cfg$dpi)
}

# --- 5. summary forest: all genes x subsets; 95% CI and prediction interval ----
sf <- pooled[Measure == "g"]
sf[, `:=`(Subset = factor(Subset, levels = c("All units", "Glomerular", "Tubulointerstitial",
                                             "Whole biopsy")),
          Gene = factor(Gene, levels = rev(cfg$genes)),
          sig  = fifelse(FDR < 0.05, "FDR < 0.05", "n.s."),
          lab  = sprintf("%.2f [%.2f, %.2f]  k=%d%s", estimate, ci_lb, ci_ub, k,
                         fifelse(is.na(I2), "", sprintf(", I²=%.0f%%", I2))))]
p <- ggplot(sf, aes(y = Gene)) +
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
ggsave(file.path(cfg$out$pool, "pooling_summary_forest.png"), p, width = 13, height = 4.2,
       dpi = cfg$dpi)

# --- 6. indication-stratified effects (Plot 4): each indication vs its own
#        study's controls, Hedges' g, DL-pooled across units (SUBSETS left out,
#        as every table and cell here is pooled) ----------------------------
pool_unc <- unc[!StudyID %in% SUBSETS]
li <- merge(pool_unc, indic[, .(StudyID, GSM, Indication)], by = c("StudyID", "GSM"))
stopifnot(nrow(li) == nrow(pool_unc))
ctrl <- li[Indication == "Control", .(m_c = mean(Expression), s_c = sd(Expression), n_c = .N),
           by = .(Gene, Batch, Compartment)]
dis  <- li[Indication != "Control", .(m_d = mean(Expression), s_d = sd(Expression), n_d = .N),
           by = .(Gene, Batch, Compartment, Indication)]
ie <- merge(dis, ctrl, by = c("Gene", "Batch", "Compartment"))[
  n_d >= pmin_n & n_c >= pmin_n & s_d > 0 & s_c > 0]
es_i <- escalc(measure = "SMD", m1i = m_d, sd1i = s_d, n1i = n_d,
               m2i = m_c, sd2i = s_c, n2i = n_c, data = ie)
ie[, `:=`(g = as.numeric(es_i$yi), g_var = as.numeric(es_i$vi))]
by_ind <- ie[, {
  if (.N >= 2) {
    f <- rma(yi = g, vi = g_var, data = data.frame(g, g_var), method = cfg$pool$method)
    list(k = .N, n_disease = sum(n_d), estimate = as.numeric(f$beta), ci_lb = f$ci.lb,
         ci_ub = f$ci.ub, p = f$pval, I2 = f$I2, units = paste(sort(Batch), collapse = ", "))
  } else with(wald_ci(g, g_var),
    list(k = 1L, n_disease = n_d, estimate = g, ci_lb = lb, ci_ub = ub,
         p = 2 * pnorm(-abs(g / sqrt(g_var))), I2 = NA_real_, units = Batch))
}, by = .(Gene, Compartment, Indication)]
by_ind[, FDR := p.adjust(p, "BH")]
write_xlsx(sig4(by_ind), file.path(cfg$out$pool, "pooling_by_indication.xlsx"))
dropped <- setdiff(unique(li[Indication != "Control", paste(Batch, Indication)]),
                   unique(ie[, paste(Batch, Indication)]))
if (length(dropped))
  msg("indication arms with < ", pmin_n, " samples or no in-study controls, left out: ",
      paste(dropped, collapse = "; "))

hm <- copy(by_ind)
hm[, `:=`(Gene = factor(Gene, levels = cfg$genes),
          row  = sprintf("%s (k=%d, n=%d)", Indication, k, n_disease),
          txt  = sprintf("%.2f%s", estimate, fifelse(FDR < 0.05, "*", "")))]
hm[, row := factor(row, levels = rev(unique(row[order(Compartment, -n_disease)])))]
lim <- 2
p <- ggplot(hm, aes(Gene, row, fill = pmax(pmin(estimate, lim), -lim))) +
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
ggsave(file.path(cfg$out$pool, "pooling_indication_heatmap.png"), p, width = 8.5,
       height = 1.4 + 0.3 * uniqueN(hm$row), dpi = cfg$dpi)

# --- 7. sensitivity to the control tissue ----------------------------------------
# Glomerular units only (the one compartment with enough units of both types),
# Living donor (k = 4, incl. GSE93798) vs Tumor nephrectomy (k = 5, incl. GSE141295
# RNA-Seq) as a moderator. ERCB holds both types; its controls are restricted to
# the living donors so each unit has one type. "Not stated" (GSE30122, DKD only)
# and "Reference RNA" (GSE66494) are not tested: each coincides with one disease.
# GSE142025 (tumor nephrectomy) and GSE111154 (Deceased donor) are whole biopsy,
# GSE35488 (Living donor) and GSE175759 (Nephrectomy) tubulointerstitial, so
# none is in the glomerular test. Per-unit effects (study_effects
# table) cover all units; the pooling and the test leave SUBSETS out.
stopifnot(all(unique(unc$StudyID) %in% names(CONTROL_TYPE)))
lc <- merge(unc, indic[, .(StudyID, GSM, Label)], by = c("StudyID", "GSM"))
lc[, ControlType := CONTROL_TYPE[StudyID]]
lc <- lc[!(Condition == "Control" & ControlType == "Living donor" &
           grepl("tumor nephrectom", Label, ignore.case = TRUE))]
cw <- arm_summaries(lc, c("Gene", "Batch", "Compartment", "ControlType"))
es_c <- effect_sizes(cw, "SMD")
cw[, `:=`(g = as.numeric(es_c$yi), g_var = as.numeric(es_c$vi))]
cw[, Pooled := !Batch %in% SUBSETS]
write_xlsx(cw[, .(Gene, Batch, Compartment, ControlType, n_Control, n_Disease,
                  g = signif(g, 4), g_var = signif(g_var, 4), Pooled)],
           file.path(cfg$out$pool, "pooling_control_type_study_effects.xlsx"))

ct_pool <- list(); ct_test <- list()
for (gene in cfg$genes) {
  d <- cw[Gene == gene & Compartment == "Glomerular" & Pooled]
  for (ct in unique(d$ControlType)) {
    x <- d[ControlType == ct]
    f <- if (nrow(x) >= 2) rma(yi = g, vi = g_var, data = x, method = cfg$pool$method)
    ci <- if (is.null(f)) wald_ci(x$g, x$g_var) else list(lb = f$ci.lb, ub = f$ci.ub)
    ct_pool[[length(ct_pool) + 1]] <- data.table(Gene = gene, ControlType = ct, k = nrow(x),
      estimate = if (is.null(f)) x$g else as.numeric(f$beta), ci_lb = ci$lb, ci_ub = ci$ub,
      I2 = if (is.null(f)) NA_real_ else f$I2)
  }
  ga <- d[ControlType %in% c("Living donor", "Tumor nephrectomy")]
  ga[, ControlType := factor(ControlType, levels = c("Living donor", "Tumor nephrectomy"))]
  fa <- rma(yi = g, vi = g_var, mods = ~ ControlType, data = ga, method = cfg$pool$method)
  ct_test[[gene]] <- data.table(Gene = gene,
    glom_k_LD = sum(ga$ControlType == "Living donor"),
    glom_k_TN = sum(ga$ControlType == "Tumor nephrectomy"),
    glom_diff_TN_minus_LD = as.numeric(fa$beta[2]), glom_diff_ci_lb = fa$ci.lb[2],
    glom_diff_ci_ub = fa$ci.ub[2], glom_QM_p = fa$QMp)
}
ct_pool <- rbindlist(ct_pool)
ct_test <- rbindlist(ct_test)
ct_test[, glom_FDR := p.adjust(glom_QM_p, "BH")]
write_xlsx(sig4(ct_pool), file.path(cfg$out$pool, "pooling_control_type_glomerular.xlsx"))
write_xlsx(sig4(ct_test), file.path(cfg$out$pool, "pooling_control_type_test.xlsx"))
msg("control-type test (glomerular, TN - LD): smallest FDR ", signif(min(ct_test$glom_FDR), 2))

cp <- copy(ct_pool)
cp[, `:=`(Gene = factor(Gene, levels = rev(cfg$genes)),
          ControlType = factor(ControlType, levels = c("Living donor", "Tumor nephrectomy",
                                                       "Not stated")))]
glab <- ct_test[, setNames(sprintf("%s\nTN − LD %.2f, p = %.2f", Gene,
                                   glom_diff_TN_minus_LD, glom_QM_p), Gene)]
p <- ggplot(cp, aes(estimate, Gene, colour = ControlType)) +
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
ggsave(file.path(cfg$out$pool, "pooling_control_type_sensitivity.png"), p, width = 9,
       height = 5.5, dpi = cfg$dpi)
write_session(cfg$out$pool, "pooling")


# #############################################################################
# PART G — MEGA-ANALYSIS WITHIN PLATFORM TIERS (plan §5)
#   "... complemented by mega-analysis within platform tiers." The samples of
#   every study on one tier (array vendor) are pooled into one model per gene:
#     nlme::lme(Expression ~ 0 + Compartment + Compartment:Disease,
#               random = pdDiag(~ Disease) | Batch, weights = varIdent(~ 1 | Batch))
#   Random study level and disease effect (the sample-level analogue of tau²),
#   residual variance per study; falls back to a random intercept if the slope
#   does not converge. A one-study tier is that study's lm(Expression ~ Disease).
#   Only Affymetrix pools studies (Illumina, Agilent, RNA-Seq = GSE142025: one
#   each; no C4BPA in the RNA-Seq tier). log2 scale;
#   a within-study z-score scale was tried and dropped (it collapses the
#   between-study disease variance to 0). Aggregated, so SUBSETS are left out
#   (this was the former sensitivity run; it is now the only analysis).
# #############################################################################
set.seed(1)
msg("Part G — mega-analysis")
mg <- copy(unc[!StudyID %in% SUBSETS])
mg[, `:=`(Tier = unname(VENDOR[vapply(COHORTS[StudyID], `[[`, "", "idtype")]),
          Disease = as.integer(Condition == "Disease"), Compartment = factor(Compartment))]
stopifnot(!anyNA(mg$Tier))
tiers <- mg[, .(units = uniqueN(Batch), studies = uniqueN(StudyID), samples = uniqueN(GSM),
                compartments = paste(sort(unique(as.character(Compartment))), collapse = ", ")),
            by = Tier]
write_xlsx(tiers, file.path(cfg$out$mega, "mega_tiers.xlsx"))
msg("platform tiers: ", paste(tiers[, sprintf("%s %d studies", Tier, studies)], collapse = ", "))

# A compartment with ONE study in the tier gets that study's own model: in a joint
# lme its random disease slope would be borrowed from the other compartments.
mega_fit <- function(d) {
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
  ctrl <- nlme::lmeControl(maxIter = 200, msMaxIter = 200, niterEM = 50, returnObject = FALSE)
  fx   <- Expression ~ 0 + Compartment + Compartment:Disease
  wts  <- nlme::varIdent(form = ~ 1 | Batch)
  model <- "lme, random intercept + disease slope"
  f <- tryCatch(nlme::lme(fx, random = list(Batch = nlme::pdDiag(~ Disease)), weights = wts,
                          data = d, method = "REML", control = ctrl),
                error = function(e) NULL)
  if (is.null(f)) {
    f <- nlme::lme(fx, random = ~ 1 | Batch, weights = wts, data = d, method = "REML", control = ctrl)
    model <- "lme, random intercept (slope did not converge)"
  }
  tt  <- summary(f)$tTable
  dis <- grep(":Disease$", rownames(tt))
  cps <- sub("^Compartment(.*):Disease$", "\\1", rownames(tt)[dis])
  vc  <- nlme::VarCorr(f)
  rbindlist(c(out, list(data.table(Compartment = cps,
             k = vapply(cps, function(cp) uniqueN(d[Compartment == cp, Batch]), 1L),
             n = vapply(cps, function(cp) nrow(d[Compartment == cp]), 1L),
             estimate = tt[dis, "Value"], SE = tt[dis, "Std.Error"], df = tt[dis, "DF"],
             p = tt[dis, "p-value"],
             tau_disease = if ("Disease" %in% rownames(vc)) as.numeric(vc["Disease", "StdDev"])
                           else NA_real_,
             model = model))))
}
run_mega <- function(d, analysis) {
  out <- list()
  for (tier in unique(d$Tier)) for (gene in cfg$genes) if (nrow(d[Tier == tier & Gene == gene]))
    out[[length(out) + 1]] <- cbind(data.table(Analysis = analysis, Tier = tier, Gene = gene),
                                    mega_fit(d[Tier == tier & Gene == gene]))
  out <- rbindlist(out)
  out[, `:=`(ci_lb = estimate - qt(0.975, df) * SE, ci_ub = estimate + qt(0.975, df) * SE)]
  out[, FDR := p.adjust(p, "BH"), by = .(Analysis, Tier, Compartment)]   # across the 5 genes
  out[]
}
primary <- sprintf("Primary (%d datasets)", uniqueN(mg$StudyID))
mega <- run_mega(mg, primary)
setcolorder(mega, c("Analysis", "Tier", "Gene", "Compartment", "k", "n", "estimate", "SE",
                    "ci_lb", "ci_ub", "df", "p", "FDR", "tau_disease", "model"))
num <- names(mega)[vapply(mega, is.numeric, TRUE) & !names(mega) %in% c("k", "n", "df")]
write_xlsx(copy(mega)[, (num) := lapply(.SD, signif, 4), .SDcols = num],
           file.path(cfg$out$mega, "mega_estimates.xlsx"))
if (any(grepl("did not converge", mega$model)))
  msg("random-slope fallback used for: ",
      paste(unique(mega[grepl("did not converge", model), paste(Analysis, Tier, Gene)]), collapse = "; "))

# Mega vs meta: the same contrast (log2FC) by two pooling routes. Glomerular and
# tubulointerstitial meta-analyses pool exactly the Affymetrix-tier units.
meta_fc <- pooled[Measure == "log2FC" & Subset != "All units",
                  .(Gene, Compartment = Subset, Tier = "Meta-analysis (all tiers)", k,
                    estimate, ci_lb, ci_ub, FDR)]
mvm <- rbind(mega[Analysis == primary, .(Gene, Compartment, Tier = paste0("Mega: ", Tier), k,
                                          estimate, ci_lb, ci_ub, FDR)], meta_fc)
write_xlsx(sig4(mvm), file.path(cfg$out$mega, "mega_vs_meta.xlsx"))
mvm_ag <- merge(mega[Analysis == primary & Tier == "Affymetrix", .(Gene, Compartment, mega = estimate)],
                meta_fc[, .(Gene, Compartment, meta = estimate)], by = c("Gene", "Compartment"))
mvm_stats <- mvm_ag[, .(n = .N, pearson_r = round(cor(mega, meta), 3),
                        same_sign = mean(sign(mega) == sign(meta)),
                        max_abs_diff = round(max(abs(mega - meta)), 3))]
write_xlsx(mvm_stats, file.path(cfg$out$mega, "mega_vs_meta_agreement.xlsx"))
msg("Affymetrix mega vs meta: r = ", mvm_stats$pearson_r, ", max |diff| = ", mvm_stats$max_abs_diff)

pd <- copy(mvm)
pd[, `:=`(Gene = factor(Gene, levels = rev(cfg$genes)),
          Compartment = factor(Compartment, levels = c("Glomerular", "Tubulointerstitial", "Whole biopsy")),
          Tier = factor(Tier, levels = c("Meta-analysis (all tiers)", "Mega: Affymetrix",
                                         "Mega: Illumina", "Mega: Agilent", "Mega: RNA-Seq")),
          sig = fifelse(FDR < 0.05, "FDR < 0.05", "n.s."))]
# group = Tier on every layer, so the significance fill does not regroup the
# points and dodge them out of line with their error bars.
p <- ggplot(pd, aes(estimate, Gene, colour = Tier, shape = Tier, group = Tier)) +
  geom_vline(xintercept = 0, colour = "grey55", linetype = 2) +
  geom_errorbar(aes(xmin = ci_lb, xmax = ci_ub), width = 0.25, orientation = "y",
                position = position_dodge(0.6)) +
  geom_point(aes(fill = sig), size = 2.6, position = position_dodge(0.6)) +
  facet_wrap(~ Compartment, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = c("Meta-analysis (all tiers)" = "#1F3B5C", "Mega: Affymetrix" = "#1B7F79",
                                 "Mega: Illumina" = "#B2182B", "Mega: Agilent" = "#E08214",
                                 "Mega: RNA-Seq" = "#7B3294"), name = NULL) +
  scale_shape_manual(values = c(23, 21, 24, 22, 25), name = NULL) +
  scale_fill_manual(values = c("FDR < 0.05" = "grey20", "n.s." = "white"), name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(title = "Mega-analysis within platform tiers vs random-effects meta-analysis",
       subtitle = paste("log2 fold change, Disease − Control, 95% CI. Meta: DL-pooled per-study log2FC;",
                        "mega: one mixed model per tier (random study level and disease effect)"),
       x = "log2 fold change", y = NULL,
       caption = paste0(sprintf("Affymetrix tier: glomerular k = %d, tubulointerstitial k = %d. ",
                                uniqueN(mg[Tier == "Affymetrix" & Compartment == "Glomerular", Batch]),
                                uniqueN(mg[Tier == "Affymetrix" & Compartment == "Tubulointerstitial", Batch])),
                        "RNA-Seq: GSE141295 (glomerular), GSE175759 (tubulointerstitial), GSE142025 (whole biopsy). ",
                        "Whole biopsy: Affymetrix = GSE111154, Illumina = GSE116626, Agilent = GSE66494.\n",
                        "Filled = FDR < 0.05.")) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 8.5),
        plot.caption = element_text(colour = "grey45", size = 7, hjust = 0),
        panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        axis.text.y = element_text(face = "bold"), legend.position = "top")
ggsave(file.path(cfg$out$mega, "mega_vs_meta.png"), p, width = cfg$mega$width,
       height = cfg$mega$height, dpi = cfg$dpi)
write_session(cfg$out$mega, "mega")


# #############################################################################
# PART H — SUMMARY FIGURES
#   1. Method agreement: Stage 4 (one WLS model, sex-adjusted) vs the pooling
#      strategy (per-study effects, DL-pooled) — the same Disease − Control log2
#      contrast by two independent routes, per study unit and per compartment.
#   2. PRISMA-style flow of dataset selection: 16 datasets in the curated
#      inventory (= plan Annex 2; the plan's "32" is not the selection count)
#      + 8 added to Raw_filtered on 2026-11-03 -> 9 excluded (2 no control arm,
#      4 duplicates, 2 no expression values, 1 too few patients) -> 15 included
#      (14 microarray + GSE142025 RNA-Seq).
# #############################################################################
set.seed(1)
msg("Part H — summary figures")

# Per study unit (single-dataset estimates, all 15 units): Stage 4's within-unit
# contrast vs the plain difference of means the pooling strategy starts from.
# Per compartment: the two aggregated estimates (both without SUBSETS).
meth <- rbind(
  merge(per_study[, .(Gene = as.character(Gene), Batch, Compartment, s4 = log2FC,
                      s4_lo = log2FC - 1.96 * SE, s4_hi = log2FC + 1.96 * SE)],
        eff[, .(Gene = as.character(Gene), Batch, pl = log2FC,
                pl_lo = log2FC - 1.96 * sqrt(log2FC_var), pl_hi = log2FC + 1.96 * sqrt(log2FC_var))],
        by = c("Gene", "Batch"))[, Level := "Per study unit"],
  merge(res[Family == "Disease", .(Gene, Compartment = sub("^Disease vs Control \\| ", "", Term),
                                   s4 = log2FC, s4_lo = log2FC - 1.96 * SE, s4_hi = log2FC + 1.96 * SE)],
        pooled[Measure == "log2FC" & Subset != "All units",
               # 4 significant figures: summary_figures.R reads the rounded pooling_estimates.xlsx
               .(Gene, Compartment = Subset, pl = signif(estimate, 4), pl_lo = signif(ci_lb, 4),
                 pl_hi = signif(ci_ub, 4))],
        by = c("Gene", "Compartment"))[, Level := "Pooled per compartment"],
  fill = TRUE)
meth[, `:=`(Gene = factor(Gene, levels = cfg$genes),
            Level = factor(Level, levels = c("Per study unit", "Pooled per compartment")))]
meth_stats <- meth[, .(n = .N, pearson_r = cor(s4, pl), slope = coef(lm(pl ~ s4))[2],
                       same_sign = mean(sign(s4) == sign(pl)), max_abs_diff = max(abs(s4 - pl))),
                   by = Level]
write_xlsx(meth, file.path(cfg$out$summary, "method_agreement.xlsx"))
write_xlsx(sig4(meth_stats), file.path(cfg$out$summary, "method_agreement_stats.xlsx"))
msg("Stage 4 vs pooling: ", meth_stats[, paste(sprintf("%s r = %.3f", Level, pearson_r), collapse = "; ")])

facet_lab <- meth_stats[, .(Level, lab = sprintf("%s\nn = %d, r = %.3f, same sign %.0f%%",
                                                 Level, n, pearson_r, 100 * same_sign))]
meth[facet_lab, Facet := i.lab, on = "Level"]
meth[, Facet := factor(Facet, levels = facet_lab$lab)]
p <- ggplot(meth, aes(s4, pl, colour = Gene, shape = Compartment)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = 2) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_errorbar(data = meth[Level == "Pooled per compartment"],
                aes(ymin = pl_lo, ymax = pl_hi), width = 0, alpha = 0.5) +
  geom_errorbar(data = meth[Level == "Pooled per compartment"],
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
ggsave(file.path(cfg$out$summary, "method_agreement.png"), p, width = 11, height = 5.2, dpi = cfg$dpi)

# PRISMA-style flow.
n_inc <- uniqueN(unc$StudyID)
n_rna <- sum(vapply(COHORTS, `[[`, "", "idtype") == "symbol")   # included RNA-Seq series
smp   <- unique(unc[!StudyID %in% SUBSETS, .(StudyID, GSM, Batch, Compartment, Condition)])
comp  <- smp[, .(n = .N, k = uniqueN(Batch)), by = Compartment][order(-n)]
cond  <- smp[, .N, by = Condition]
excl  <- cfg$summary$excluded
stopifnot(n_inc + length(excl) == cfg$summary$n_listed)
ex_n  <- table(factor(excl, levels = unique(excl)))
write_xlsx(data.table(Step  = c("Listed (inventory + added 2026-11-03)", paste("Excluded:", names(ex_n)), "Included"),
                      Count = c(cfg$summary$n_listed, as.integer(ex_n), n_inc)),
           file.path(cfg$out$summary, "prisma_counts.xlsx"))
boxes <- data.table(
  x = c(0, 0, 1.3), y = c(3, 1.45, 2.45), w = c(0.95, 0.95, 1.25), h = c(0.55, 1.4, 0.95),
  fill = c("main", "inc", "exc"),
  label = c(
    sprintf("GEO datasets: curated inventory (16)\n+ added 2026-11-03 (8),  n = %d", cfg$summary$n_listed),
    sprintf(paste0("Included in the analysis: n = %d series\n(%d microarray + %d RNA-Seq)\n",
                   "Pooled analyses: %d series, %d units, %d samples\n",
                   "(GSE30122 split by compartment)\n",
                   "%s\n",
                   "Control %d / Disease %d;  %d disease indications"),
            n_inc, n_inc - n_rna, n_rna, uniqueN(smp$StudyID), uniqueN(smp$Batch),
            uniqueN(smp$GSM),
            paste(sprintf("%s: %d samples, k = %d", comp$Compartment, comp$n, comp$k), collapse = "\n"),
            cond[Condition == "Control", N], cond[Condition == "Disease", N],
            uniqueN(indic[!StudyID %in% SUBSETS & Indication != "Control", Indication])),
    paste0(sprintf("Excluded: n = %d\n", length(excl)),
           paste(sprintf("%s: %s", names(ex_n),
                         vapply(names(ex_n), function(r) paste(names(excl)[excl == r], collapse = ", "), "")),
                 collapse = "\n"))))
# Listed -> included down the main column; the exclusion box branches off it.
flow <- data.table(x = c(0, 0), xend = c(0, 1.3 - 0.625),
                   y = c(3 - 0.275, 2.45), yend = c(1.45 + 0.7, 2.45))
p <- ggplot() +
  geom_rect(data = boxes, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2,
                              ymax = y + h / 2, fill = fill), colour = "grey30", linewidth = 0.4) +
  geom_text(data = boxes, aes(x, y, label = label), size = 3, lineheight = 1.05) +
  geom_segment(data = flow, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(2.2, "mm"), type = "closed"), linewidth = 0.4) +
  scale_fill_manual(values = c(main = "#E8EEF5", inc = "#DCEFEA", exc = "#F6E3E3"), guide = "none") +
  coord_cartesian(xlim = c(-0.55, 2.0), ylim = c(0.65, 3.4), expand = FALSE) +
  labs(title = "Dataset selection (PRISMA-style)") +
  theme_void(base_size = 10) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.background = element_rect(fill = "white", colour = NA))
ggsave(file.path(cfg$out$summary, "prisma_flow.png"), p, width = 10.5, height = 4.8, dpi = cfg$dpi)
write_session(cfg$out$summary, "summary")


# #############################################################################
# PART I — ML INPUT WORKBOOKS (from create_ML_train_data.R; the RNA-Seq series
#   without controls inline from rnaseq_panel.R)
#   One row per GSM, before cross-study batch correction: the five panel genes,
#   metadata, Sex (Stage 4: reported, else inferred; M / F / NA), Platform and
#   Disease_Label (Control = 0, Disease = 1). Aggregated, so SUBSETS are left
#   out as series; their samples remain through GSE30122. Microarray + RNA-Seq
#   in one table (by instruction, 2026-11-03): GSE142025 comes through Stage 2/3
#   like every cohort; GSE137570-GPL23227 and GSE219185 (no controls) enter the
#   ML inputs only. No scaling, overlap exclusion or train/test split:
#   preprocessing belongs inside the CV folds.
#   all_bn.xlsx: the same rows, the Stage 3 model refitted gene by gene over all
#   ML rows. That model protected Disease, so the labels shaped those values
#   (label leakage).
# #############################################################################
set.seed(1)
msg("Part I — ML input workbooks")

# --- RNA-Seq series WITHOUT controls (rnaseq_panel.R): ML inputs only ---------
# cond_field = characteristic holding the diagnosis; indication maps its value.
RNA_NOCTRL <- list(
  list(id = "GSE137570-GPL23227", gse = "GSE137570", idtype = "symbol",
       cond_field = "tissue", compartment = "Whole biopsy",
       indication = c(".*" = "CKD (unspecified)")),
  list(id = "GSE219185",          gse = "GSE219185", idtype = "ensembl",
       cond_field = "disease", compartment = "Tubulointerstitial",
       indication = c("focal segmental" = "FSGS", "minimal change" = "MCD")))
ML_RNASEQ <- c(names(Filter(function(co) co$idtype == "symbol", COHORTS)), vapply(RNA_NOCTRL, `[[`, "", "id"))   # Platform = "RNA-Seq"

read_pheno_ch <- function(path) {
  lines <- grep("^!Sample_(geo_accession|characteristics_ch1)\t", readLines(path, warn = FALSE),
                value = TRUE)
  rows  <- lapply(strsplit(lines, "\t", fixed = TRUE), function(p) gsub('^"|"$', "", p[-1]))
  gsm   <- rows[[which(startsWith(lines, "!Sample_geo_accession"))]]
  ch    <- rows[startsWith(lines, "!Sample_characteristics_ch1")]
  key   <- vapply(ch, function(v) tolower(sub(":.*", "", v[1])), "")
  list(gsm = gsm, ch = setNames(lapply(ch, function(v) trimws(sub("^[^:]*:", "", v))), key))
}
# Same rule as Stage 4: score = mean z(Y genes) - z(XIST), 2-means split,
# accepted only when the clusters are well separated.
infer_sex_rows <- function(ex) {
  y <- intersect(cfg$s4$y_genes, rownames(ex)); has_x <- cfg$s4$x_gene %in% rownames(ex)
  if (!length(y) && !has_x) return(rep(NA_character_, ncol(ex)))
  z <- t(scale(t(ex[c(y, if (has_x) cfg$s4$x_gene), , drop = FALSE])))
  score <- (if (length(y)) colMeans(z[y, , drop = FALSE]) else 0) -
           (if (has_x) z[cfg$s4$x_gene, ] else 0)
  km  <- kmeans(score, centers = range(score))
  sep <- abs(diff(km$centers[, 1])) /
         sqrt(mean(tapply(score, km$cluster, function(s) if (length(s) > 1) var(s) else 0)))
  if (!is.finite(sep) || sep < cfg$s4$min_sex_sep) return(rep(NA_character_, ncol(ex)))
  ifelse(km$cluster == which.max(km$centers[, 1]), "M", "F")
}
rna_long <- list(); rna_smp <- list(); rna_qc <- list()
for (co in RNA_NOCTRL) {
  dir <- file.path(cfg$raw_dir, co$gse)
  ph  <- read_pheno_ch(file.path(dir, paste0(co$id, "_series_matrix.txt")))
  lab <- ph$ch[[co$cond_field]]; stopifnot(!is.null(lab))
  ex  <- fread(file.path(dir, paste0(co$id, "_expression.txt")))
  ids <- ex[[1]]; ex <- as.matrix(ex[, -1])
  if (co$idtype == "ensembl")   # Ensembl -> symbol, only for the genes used here
    ids <- unname(suppressMessages(AnnotationDbi::mapIds(
      getExportedValue("org.Hs.eg.db", "org.Hs.eg.db"), keys = ids,
      column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")))
  keep <- ids %in% c(cfg$genes, SEX_GENES)
  stopifnot(!anyDuplicated(ids[keep]))   # else a MaxMean collapse would be needed
  ex <- ex[keep, , drop = FALSE]; rownames(ex) <- ids[keep]
  stopifnot(setequal(colnames(ex), ph$gsm))
  ex <- ex[, ph$gsm, drop = FALSE]
  rep_sex <- ph$ch[["sex"]]
  rep_sex <- if (is.null(rep_sex)) rep(NA_character_, length(ph$gsm)) else
    c(male = "M", female = "F")[tolower(rep_sex)]
  sex <- ifelse(is.na(rep_sex), infer_sex_rows(ex), rep_sex)
  ind <- rep(NA_character_, length(lab))
  for (pat in names(co$indication))
    ind[is.na(ind) & grepl(pat, lab, ignore.case = TRUE)] <- co$indication[[pat]]
  stopifnot(!anyNA(ind))
  rna_smp[[co$id]] <- data.table(StudyID = co$id, GSM = ph$gsm, Compartment = co$compartment,
    Condition = "Disease", Indication = ind, Label = lab, Sex = unname(sex))
  g <- intersect(cfg$genes, rownames(ex))
  rna_long[[co$id]] <- data.table(GSM = rep(ph$gsm, each = length(g)), StudyID = co$id,
    Condition = "Disease", Compartment = co$compartment, Batch = co$id,
    Gene = g, Expression = as.numeric(ex[g, ]))
  rna_qc[[co$id]] <- data.table(Series = co$id, n = length(ph$gsm),
    Use = "ML inputs only (no control arm)",
    Indications = paste(names(table(ind)), table(ind), collapse = ", "),
    SexSource = if (all(!is.na(rep_sex))) "reported" else "inferred",
    PanelGenes = paste(g, collapse = ", "),
    MissingPanelGenes = paste(setdiff(cfg$genes, g), collapse = ", "))
}
rna_long <- rbindlist(rna_long); rna_smp <- rbindlist(rna_smp); rna_qc <- rbindlist(rna_qc)
stopifnot(nrow(rna_long) > 0, !anyNA(rna_long$Expression))
# The same three tables rnaseq_panel.R writes.
write_xlsx(rna_long, file.path(cfg$out$rnaseq, "panel_rnaseq_long.xlsx"))
write_xlsx(rna_smp,  file.path(cfg$out$rnaseq, "rnaseq_samples.xlsx"))
write_xlsx(rna_qc,   file.path(cfg$out$rnaseq, "rnaseq_qc.xlsx"))

ml_unc <- rbind(unc[!StudyID %in% SUBSETS], rna_long)
ml_ind <- rbind(indic[!StudyID %in% SUBSETS],
                rna_smp[, .(StudyID, GSM, Compartment, Condition, Indication, Label)])
ml_sex <- rbind(sex_tab[!StudyID %in% SUBSETS, .(StudyID, GSM, Sex)], rna_smp[, .(StudyID, GSM, Sex)])
stopifnot(
  setequal(unique(ml_unc$Gene), cfg$genes),
  all(is.finite(ml_unc$Expression)),
  all(ml_unc$Condition %in% c("Control", "Disease")),
  !anyDuplicated(ml_unc[, .(StudyID, GSM, Gene)]),
  !anyDuplicated(ml_ind[, .(StudyID, GSM)]),
  !anyDuplicated(ml_sex[, .(StudyID, GSM)]),
  all(is.na(ml_sex$Sex) | ml_sex$Sex %in% c("M", "F"))
)

training_data <- dcast(ml_unc, GSM + StudyID + Batch + Compartment + Condition ~ Gene,
                       value.var = "Expression")
training_data <- merge(training_data, ml_ind[, .(StudyID, GSM, Indication)],
                       by = c("StudyID", "GSM"), all.x = TRUE, sort = FALSE)
stopifnot(nrow(training_data[!ml_sex, on = .(StudyID, GSM)]) == 0L)
training_data <- merge(training_data, ml_sex[, .(StudyID, GSM, Sex)],
                       by = c("StudyID", "GSM"), all.x = TRUE, sort = FALSE)
training_data[, Disease_Label := as.integer(Condition == "Disease")]
training_data[, Platform := fifelse(StudyID %in% ML_RNASEQ, "RNA-Seq", "Microarray")]
ml_cols <- c("GSM", "StudyID", "Batch", "Compartment", "Condition",
             "Disease_Label", "Indication", "Sex", "Platform", cfg$genes)
setcolorder(training_data, ml_cols)
setorder(training_data, StudyID, GSM)
# A gene is NA only where the study's deposit lacks it (C4BPA: GSE142025,
# GSE137570); every value present is finite.
present <- unique(ml_unc[, .(StudyID, Gene)])
na_ok <- function(g) !nrow(present[Gene == g][training_data[is.na(get(g))], on = "StudyID", nomatch = 0L])
stopifnot(
  nrow(training_data) == 1322L,                 # 931 microarray + 36 + 24 + 87 + 24 + 220 RNA-Seq
  !anyDuplicated(training_data$GSM),
  uniqueN(training_data$StudyID) == 19L,
  sum(training_data$Disease_Label == 1L) == 1077L,
  sum(training_data$Disease_Label == 0L) == 245L,
  all(vapply(cfg$genes, na_ok, TRUE)),
  all(vapply(training_data[, cfg$genes, with = FALSE],
             function(x) all(is.finite(x[!is.na(x)])), TRUE))
)

ml_readme <- data.table(
  Item = c("Rows", "Predictors", "Sex", "Sex provenance", "Target", "Metadata", "Expression scale",
           "Batch correction", "Excluded series", "Sample independence",
           "Cross-study overlap", "Validation", "Input provenance",
           "Platform", "Missing values", "Studies without controls"),
  Description = c(
    "1322 unique GEO GSM accessions: 1077 disease and 245 control records (931 microarray, 391 RNA-Seq).",
    paste0(paste(cfg$genes, collapse = ", "), "; Sex is an optional categorical predictor."),
    "Sex: M = Male; F = Female; missing values remain NA. This is recorded/inferred sex, not gender identity.",
    "Reported sex takes priority; otherwise inferred by the Stage4 XIST/Y-marker procedure, not the five RenoRisk genes. Per-sample provenance (reported/inferred/unknown) is in Results/Stage4/sex_by_sample.xlsx (GSE137570/GSE219185: Results/RNASeq/rnaseq_qc.xlsx), not in this workbook.",
    "Disease_Label: 0 = Control; 1 = Disease. Condition preserves the source label.",
    "GSM, StudyID, Batch, Compartment, Indication and Platform describe samples; they are not model predictors.",
    "Stage23 uncorrected values: deposited log2/relative expression; GSE66494 uses log2(x + 1). Gene-centred cohorts retain their relative scale. RNA-Seq rows are the submitters' normalised log2 tables (GSE219185 voom log2-CPM; GSE141295 and GSE175759 DESeq2 VST computed here from the deposited counts by rnaseq_counts_vst.R; GSE142025 and GSE137570 unstated): a different scale and dynamic range from the arrays.",
    "No cross-study limma correction or additional ML scaling is applied to this workbook.",
    paste("GSE30528 and GSE30529 (all samples in GSE30122) and GSE35487 (same arrays as GSE35488) and GSE37460 (arrays already in GSE104948/GSE99339) are not used."),
    "Unique GSM accessions are not necessarily independent biological samples or patients. Documented reanalyses remain included.",
    "The expression audit flags GSE108109/GSE200828 as strongly overlapping and GSE121233/GSE129973 as additional candidates. See Results/SampleOverlap/README.md. GSE219185 (NEPTUNE tubulointerstitium) states it shares samples with NEPTUNE GEO sets including GSE108113, whose glomerular part is GSE108109: same patients, different compartment.",
    "Use study/patient-aware validation (leave-one-study-out). Keep linked studies together (e.g. GSE219185 with GSE108109); fit learned preprocessing only within training folds. Sex is coded M / F (categorical). No fixed train/test split is assigned here.",
    "Expression/indications: Results/Stage23/panel_uncorrected_long.xlsx and sample_indication.xlsx. Sex: Results/Stage4/sex_by_sample.xlsx. GSE137570/GSE219185: Results/RNASeq/ (rnaseq_panel.R). Submitter-processed values from Data/Raw_filtered only.",
    "Microarray (931 rows) or RNA-Seq (391 rows: GSE142025 36, GSE141295 24, GSE175759 87, GSE137570-GPL23227 24, GSE219185 220). Platform is confounded with study: use it to stratify, not as a predictor.",
    "C4BPA is NA for GSE142025 and GSE137570 (60 rows): not in those deposits. Impute inside training folds or drop C4BPA; most scikit-learn estimators reject NaN.",
    "GSE137570 (CKD) and GSE219185 (FSGS/MCD) are all Disease. Their study identity predicts the label, so only grouped-by-study CV is valid; random-split CV is optimistic."
  ))

# Write, then read back: the workbook must preserve feature values and labels.
# Values are compared to 1e-12, not identical(): an in-memory double can lose its
# last bit on the first xlsx round trip (relative error ~5e-16).
ml_col_types <- c(rep("text", 5L), "numeric", "text", "text", "text", rep("numeric", 5L))
write_ml <- function(d, readme, file, dir = cfg$out$ml) {
  path <- file.path(dir, file)
  write_xlsx(list(training_data = d, README = readme), path)
  back <- as.data.table(readxl::read_excel(path, sheet = "training_data", col_types = ml_col_types))
  stopifnot(identical(back$GSM, d$GSM),
            isTRUE(all.equal(as.matrix(back[, cfg$genes, with = FALSE]),
                             as.matrix(d[, cfg$genes, with = FALSE]), tolerance = 1e-12)),
            identical(as.integer(back$Disease_Label), d$Disease_Label),
            identical(back$Sex, d$Sex))
  msg("ML input: ", path)
}
write_ml(training_data, ml_readme, "all.xlsx", cfg$out$ml_input)

# Same rows and metadata, batch-corrected. The Stage 3 model (batch = StudyID;
# design = intercept + Disease), refitted
# over ALL ML rows so the RNA-Seq series without controls (not in Stage 3) are
# corrected too. Gene by gene on the rows that have the gene. A disease-only
# study's offset rests on the assumption that its disease effect equals the
# pooled one.
batch_normalized <- copy(training_data)
for (g in cfg$genes) {
  ok <- !is.na(batch_normalized[[g]])
  d  <- batch_normalized[ok]
  b  <- factor(d$StudyID)
  des <- model.matrix(~ Disease_Label, data = d)
  X.batch <- model.matrix(~ b, contrasts.arg = list(b = contr.sum(levels(b))))[, -1, drop = FALSE]
  stopifnot("design + batch is rank-deficient" =
              qr(cbind(des, X.batch))$rank == ncol(des) + ncol(X.batch))
  set(batch_normalized, which(ok), g,
      as.numeric(limma::removeBatchEffect(matrix(d[[g]], nrow = 1), batch = b, design = des)))
}
stopifnot(
  identical(batch_normalized$GSM, training_data$GSM),
  identical(batch_normalized$Disease_Label, training_data$Disease_Label),
  all(vapply(cfg$genes, function(g)
    identical(is.na(batch_normalized[[g]]), is.na(training_data[[g]])), TRUE))
)

ml_readme_bn <- copy(ml_readme)
ml_readme_bn[Item == "Expression scale", Description :=
  "limma::removeBatchEffect with the Stage23 model (plan: ~ StudyID + Condition; batch = StudyID, 19 studies; design = intercept + Disease), refitted over all 1322 rows gene by gene, so values differ from Results/Stage23/panel_harmonised_long.xlsx (fitted on 1078 rows). Offsets of the two disease-only RNA-Seq series (GSE137570, GSE219185) assume their disease effect equals the pooled one. Cross-compartment level differences are removed except within GSE30122 (one study, both compartments). Gene-centred cohorts make the common scale relative, not absolute."]
ml_readme_bn[Item == "Batch correction", Description :=
  "Corrected for batch location only, not spread: residual per-gene SD still differs up to ~6.7x between studies. The correction design protected a Disease effect, so Disease_Label influenced these values."]
ml_readme_bn[Item == "Validation", Description :=
  "LABEL LEAKAGE: the correction was fitted on all samples with their disease labels, including any later test fold. Cross-validated performance on these values is optimistic. Use Data/ML_input/all.xlsx with in-fold preprocessing for unbiased validation; use this workbook for visualisation or comparison. Sex is coded M / F (categorical)."]
ml_readme_bn[Item == "Input provenance", Description :=
  "Expression: Data/ML_input/all.xlsx values, corrected in create_ML_train_data.R section 6. Submitter-processed values from Data/Raw_filtered only."]
write_ml(batch_normalized, ml_readme_bn, "all_bn.xlsx")
msg("ML input: ", paste(training_data[, .N, by = Condition][, sprintf("%s %d", Condition, N)],
                        collapse = ", "), "; sex ",
    paste(training_data[, .N, by = Sex][, sprintf("%s %d", fifelse(is.na(Sex), "NA", Sex), N)],
          collapse = ", "))
write_session(cfg$out$ml, "ML_training")

msg("all parts done in ", round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1), " s")
