#!/usr/bin/env Rscript
# =============================================================================
# MetaAnalysis — Stage 2 + 3 (visualisation route): comparable 5-gene boxplots
#
# SOURCE: CKD_Biomarker_Pooled_Analysis_Plan_Final.pdf — "PROJECT ANALYSIS PLAN &
#   TECHNICAL MASTER PLAN", RenoRisk Validation Across Multiple CKD Indications,
#   Version 2.0 (Final Analysis Plan), July 2026. Project Lead: Amela Jusic;
#   Execution Lead: Stojanko Tudjarski. Quotations below are verbatim.
#
# PLAN REQUIREMENTS IMPLEMENTED HERE
#   §3 Master Analytical Pipeline & Integration Workflow:
#     "To evaluate C4BPA, COL13A1, HABP2, CFH, and CFI across diverse public
#     datasets, the raw count/intensity data must undergo uniform harmonization
#     and batch correction to remove technical artifacts while preserving true
#     biological signal."
#   §3 Stage 1 — Raw Data Processing & Uniform Quantification:
#     "For RNA-Seq cohorts, align raw reads to GRCh38 and quantify transcripts.
#     For Affymetrix microarrays, perform RMA background correction and
#     quantiles normalization."
#   §3 Stage 2 — Target Gene Mapping:
#     "Specifically extract probe sets/transcripts corresponding to C4BPB,
#     COL13A1, HABP2, CFH, and CFI. For microarrays with multiple probes per gene,
#     collapse using WGCNA::collapseRows (selecting maximum mean expression)."
#   §3 Stage 3 — Empirical Bayes Batch Correction:
#     "Apply ComBat-seq for raw integer RNA-Seq counts and limma::removeBatchEffect
#     / sva::ComBat for normalized continuous expression, modeling ~ StudyID +
#     Condition."
#     (C4BPB -> C4BPA: replaced at the request of the Project Lead.)
#   §1 Project Roles (Execution Lead responsibilities):
#     "Processing and normalization of 32 GEO datasets (Microarray, RNA-Seq,
#     snRNA-Seq)"; "Cross-cohort harmonization using ComBat-seq and limma";
#     "Code generation, pipeline reproducibility, and figure preparation".
#   §2 Target RenoRisk Biomarker Panel Specification:
#     "The RenoRisk panel consists of 5 specific target biomarkers categorized into
#     primary diagnostic candidates and additional complement system regulators"
#     Main RenoRisk: C4BPA (C4b-Binding Protein Alpha), COL13A1 (Collagen Type XIII
#     Alpha 1 Chain), HABP2 (Hyaluronan Binding Protein 2); Additional RenoRisk:
#     CFH (Complement Factor H), CFI (Complement Factor I).
#   §4 Cohort Stratification & GEO Dataset Inventory:
#     "The 32 datasets identified in the GEO repository search have been mapped to
#     evaluate RenoRisk target performance across diverse clinical indications:
#     see table in annex 1."
#   §7 Plot 1 — PCA / UMAP Plots (Pre- vs. Post-Batch Correction):
#     "Demonstrates batch effect removal across GEO cohorts."
#
# HOW THIS SCRIPT ADDRESSES THEM (and where it deviates, with the reason)
#   * Stage 1 is NOT run: it was retired (Scripts/stage1_quantification_old.R).
#     The input is the submitter-processed values of the 16 microarray cohorts
#     in Data/Raw_filtered (five different processing pipelines) plus one RNA-Seq
#     cohort, GSE142025 (by instruction, 2026-11-03: RNA-Seq in the same bag; its
#     series matrix is metadata-only, values come from its submitter-normalised
#     log2 <id>_expression.txt) — a standing caveat on every figure. Those values
#     are normalised log2, not counts, so ComBat-seq does not apply.
#   * Datasets: the curated dataset inventory lists 16 GEO series (the same 16
#     accessions as plan Annex 2, "List of the study eligible for meta-analysis");
#     of the 3 RNA-Seq series only GSE142025 has controls, so 14 are used here, plus GSE93798, GSE35488, GSE111154 (added 2026-11-03) = 17.
#     GSE137570 / GSE219185 (no controls) enter the ML inputs only
#     (Scripts/by-claude/rnaseq_panel.R). The plan's "32" is not the selection count.
#   * Stage 2: probes are mapped per platform (custom-CDF Entrez IDs, or the
#     hta20transcriptcluster / hgu133a2 / illuminaHumanv4 / hgug4112a annotation
#     packages) and multi-probe genes are collapsed with
#     WGCNA::collapseRows(method = "MaxMean"), exactly as the plan specifies.
#     GSE66494 (Agilent) is deposited on a linear scale and is log2-transformed.
#   * Stage 3: limma::removeBatchEffect (one of the two methods the plan allows),
#     ONE model over the aggregated samples (GSE30528/GSE30529, subsets of
#     GSE30122, left out), modelled exactly as the plan says, "~ StudyID +
#     Condition" (by instruction, 2026-11-04; it replaced batch = study x
#     compartment unit with a Disease effect per compartment): batch = StudyID
#     (17), design = intercept + Disease. Compartment is nested in study for every
#     cohort but GSE30122, so cross-compartment level differences are removed
#     except inside GSE30122 (one batch, both compartments): compare within a
#     compartment only. One Disease effect is protected for all compartments.
#     ComBat-seq is applied to the two count-deposited RNA-Seq series
#     (GSE141295, GSE175759) in rnaseq_counts_vst.R, before their VST. removeBatchEffect corrects location, not spread
#     (residual SD still differs up to ~6.7x; residual_sd_ratio.xlsx).
#   * Plot 1: PCA and UMAP of the 5 panel genes before vs after correction, with
#     the batch share of PC1/PC2 (pca_batch_r2.xlsx) and the UMAP batch
#     silhouette (umap_batch_silhouette.xlsx) quantifying the removal.
#   * For plan §7 Plot 4 (indication heatmap, built in pooling_strategy.R) this
#     script also records each sample's specific diagnosis (sample_indication.xlsx).
#
# Input : the curated series matrices under Data/Raw_filtered/
# Output: Results/Stage23/ — tables (.xlsx) and figures (panel boxplots, PCA, UMAP)
#
# WHY THIS SCRIPT EXISTS
#   The 17 curated cohorts (16 arrays + 1 RNA-Seq) are NOT mutually comparable as deposited:
#     * 7 platforms (6 array + RNA-Seq), 5 ID namespaces (incl. gene symbols)
#     * 4 incompatible value scales — log2-absolute, per-gene-centred, linear, RNA-Seq normalised log2
#   So a shared Y axis requires: symbol mapping -> scale repair -> batch removal.
#
# HARMONISATION STRATEGY (Option B of the design discussion)
#   limma::removeBatchEffect(x, batch = StudyID, design = model.matrix(~Condition))
#   removes the additive per-gene study offset while PRESERVING the Control vs
#   Disease contrast. Unlike a plain per-cohort z-score it does not flatten real
#   biology in cohorts whose case/control ratio is skewed.
#
# COHORT SET AND COMPARTMENT HANDLING
#   1. All 16 microarray Raw_filtered datasets are read. GSE30122 fully contains GSE30528
#      (22/22 GSMs) and GSE30529 (22/22): those two (cfg$subsets) are kept in
#      every single-dataset table (cohort_qc, sample_indication, and the
#      uncorrected hand-off panel_uncorrected_long, from which per-study effects
#      are computed) but left out of everything aggregated across datasets —
#      the removeBatchEffect fit, the corrected values and every figure — so no
#      sample is counted twice (967 samples, 16 batch units).
#   3. A panel gene missing from a deposit (C4BPA in GSE142025) stays in the panel:
#      that cohort is NA for it, the correction is fitted gene by gene on the
#      samples that have the gene, and PCA / UMAP use the genes on every cohort.
#   2. All aggregated samples are corrected in ONE removeBatchEffect model. Compartment is
#      perfectly nested inside study for most cohorts, so compartment main
#      effects are not estimable alongside batch; they are absorbed into the
#      batch terms (batch = StudyID, design ~ Condition, as the plan says).
#      GSE30122 carries both compartments in one batch, so its compartment
#      offset stays. Consequence: cross-compartment LEVEL differences are removed
#      elsewhere — compare Control vs Disease within a compartment, not compartments.
#
# NOTE ON INPUT PROVENANCE
#   These are submitter-processed values (GeneSpring, Expression Console, a
#   "modified Affy package", Illumina quantile, Agilent FeatureExtraction — five
#   different pipelines). Stage 1 (uniform re-quantification from raw data,
#   Scripts/stage1_quantification_old.R) was meant to bypass them but is
#   retired, so every figure and table carries this limitation and must say so.
# =============================================================================

cfg <- list(
  project_dir = path.expand("~/Striki/Quant/quant.2/MetaAnalysis"),
  genes       = c("C4BPA", "COL13A1", "HABP2", "CFH", "CFI"),
# GSE30528/GSE30529 (subsets of GSE30122) were dropped and renamed -dup on 2026-11-03 (by instruction); the guard stays for any future contained dataset.
  subsets     = character(0),   # datasets fully contained in another: single-dataset use only
  min_group_n = 3,          # drop a stratum arm with fewer samples than this
  gridlines   = "vertical", # "vertical" | "horizontal" | "both" | "none"
                            #   vertical   = grey rules on the boundaries BETWEEN
                            #                groups, so each Control+Disease pair
                            #                sits in its own lane (geom_vline)
                            #   horizontal = the y gridlines from theme_bw()
  width       = 11,
  height      = 4.875,  # by-compartment figure
  pc_width    = 13,     # per-cohort figure
  pc_height   = 22,
  pca_width   = 12,     # PCA and UMAP before/after figures
  pca_height  = 5.5,
  umap_neighbors = 15,  # uwot defaults
  umap_min_dist  = 0.1,
  dpi         = 300
)
cfg$raw_dir <- file.path(cfg$project_dir, "Data", "Raw_filtered")
cfg$results_dir <- file.path(cfg$project_dir, "Results")   # everything the pipeline generates
cfg$out_dir <- file.path(cfg$results_dir, "Stage23")      # tables + provenance
cfg$fig_dir <- cfg$out_dir                                  # figures sit with the tables

suppressPackageStartupMessages({
  library(data.table); library(limma); library(ggplot2); library(AnnotationDbi)
  library(WGCNA); library(writexl)
})
set.seed(1)
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg$fig_dir, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

# Entrez IDs — the ERCB/Entrez custom CDFs address probes as "<entrez>_at"
ENTREZ <- c(C4BPA = "722", COL13A1 = "1305", HABP2 = "3026", CFH = "3075", CFI = "3426")

# ---------------------------------------------------------------------------
# 1. Cohort registry — curated from the series-matrix metadata.
#    `control` is a regex applied to the field named in `cond_field`; samples
#    matching it are Control, all others Disease. `cond_field` "diagnosis_blank"
#    is the ERCB convention where living donors carry an EMPTY diagnosis value;
#    ERCB's tumor nephrectomies are Control too (by instruction, 2026-10-09), as
#    they are in GSE99339 / GSE129973 / GSE96804.
# ---------------------------------------------------------------------------
COHORTS <- list(
  list(id = "GSE104948-GPL22945", dir = "GSE104948", idtype = "entrez",
       compartment = "Glomerular",        cond_field = "diagnosis_blank", control = "tumor nephrectom"),
  list(id = "GSE104954-GPL22945", dir = "GSE104954", idtype = "entrez",
       compartment = "Tubulointerstitial", cond_field = "diagnosis_blank", control = "tumor nephrectom"),
  list(id = "GSE108109",          dir = "GSE108109", idtype = "entrez",
       compartment = "Glomerular",        cond_field = "source",  control = "living donor"),
  list(id = "GSE200828",          dir = "GSE200828", idtype = "entrez",
       compartment = "Glomerular",        cond_field = "source",  control = "living donor"),
  list(id = "GSE99339-GPL19184",  dir = "GSE99339",  idtype = "entrez",
       compartment = "Glomerular",        cond_field = "description", control = "tumor nephrectomy"),
  # Added 2026-11-03 (by instruction): GSE93798 (IgAN glomeruli, custom CDF v19 as
  # GSE104948), GSE35488 (IgAN tubulointerstitium, Entrez CDF v10: bare Entrez IDs),
  # GSE111154 (post-mortem cortex, early DN vs non-diabetic donors, HTA 2.0).
  list(id = "GSE93798",           dir = "GSE93798",  idtype = "entrez",
       compartment = "Glomerular",        cond_field = "source",  control = "control"),
  list(id = "GSE35488",           dir = "GSE35488",  idtype = "entrez",
       compartment = "Tubulointerstitial", cond_field = "source", control = "living donor"),
  list(id = "GSE111154",          dir = "GSE111154", idtype = "hta20transcriptcluster.db",
       compartment = "Whole biopsy",      cond_field = "source",  control = "control"),
  list(id = "GSE121233",         dir = "GSE121233", idtype = "hta20transcriptcluster.db",
       compartment = "Glomerular",        cond_field = "disease state", control = "^control"),
  list(id = "GSE129973",          dir = "GSE129973", idtype = "hta20transcriptcluster.db",
       compartment = "Glomerular",        cond_field = "source",  control = "tumor nephrectom"),
  list(id = "GSE96804",           dir = "GSE96804",  idtype = "hta20transcriptcluster.db",
       compartment = "Glomerular",        cond_field = "source",  control = "tumor nephrectom"),
  list(id = "GSE30122",           dir = "GSE30122",  idtype = "hgu133a2.db",
       compartment = "split",             cond_field = "source",  control = "control kidney"),
  list(id = "GSE116626",          dir = "GSE116626", idtype = "illuminaHumanv4.db",
       compartment = "Whole biopsy",      cond_field = "subject status", control = "living donor"),
  list(id = "GSE66494",           dir = "GSE66494",  idtype = "hgug4112a.db",
       compartment = "Whole biopsy",      cond_field = "sample type", control = "normal kidney",
       log2 = TRUE),  # Agilent FeatureExtraction output is LINEAR, not log
  # RNA-Seq (by instruction, 2026-11-03: same bag as the microarrays). The series
  # matrix holds labels only; values are the submitter's normalised log2 table in
  # `expr_file`, rows addressed by gene symbol. C4BPA is not in this deposit.
  list(id = "GSE142025",          dir = "GSE142025", idtype = "symbol",
       expr_file = "GSE142025_expression.txt",
       compartment = "Whole biopsy",      cond_field = "group", control = "^control$"),
  # RNA-Seq deposited as raw counts only (added 2026-11-03, by instruction): values
  # are DESeq2 VST made by Scripts/by-claude/rnaseq_counts_vst.R from GEO <id>_RAW.tar.
  list(id = "GSE141295",          dir = "GSE141295", idtype = "symbol",
       expr_file = "GSE141295_expression.txt",
       compartment = "Glomerular",        cond_field = "tissue1", control = "normal cortex"),
  list(id = "GSE175759",          dir = "GSE175759", idtype = "symbol",
       expr_file = "GSE175759_expression.txt",
       compartment = "Tubulointerstitial", cond_field = "diagnosis", control = "^control$")
)

# ---------------------------------------------------------------------------
# 2. Series-matrix reader.  Metadata only + the expression table; no GPL fetch.
# ---------------------------------------------------------------------------
read_sm <- function(path) {
  lines <- readLines(path, warn = FALSE)
  meta  <- grep("^!Sample_", lines, value = TRUE)
  parts <- lapply(meta, function(ln) {
    p <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    list(key = sub("^!Sample_", "", p[1]), vals = gsub('^"|"$', "", p[-1]))
  })
  keys  <- make.unique(vapply(parts, `[[`, "", "key"))
  pdata <- as.data.frame(setNames(lapply(parts, `[[`, "vals"), keys),
                         stringsAsFactors = FALSE)
  rownames(pdata) <- pdata$geo_accession

  b <- which(trimws(lines) == "!series_matrix_table_begin") + 1
  e <- which(trimws(lines) == "!series_matrix_table_end")   - 1
  expr <- as.matrix(read.table(text = paste(lines[b:e], collapse = "\n"),
                               header = TRUE, sep = "\t", row.names = 1,
                               check.names = FALSE, quote = '"', comment.char = ""))
  list(pheno = pdata, expr = expr)
}

# Pull the characteristics/source/description column holding a given field.
get_field <- function(pdata, field) {
  if (field == "source")      return(tolower(pdata$source_name_ch1))
  if (field == "description") return(tolower(pdata$description))
  cols <- grep("^characteristics_ch1", names(pdata), value = TRUE)
  hit  <- cols[vapply(cols, function(c)
            any(grepl(paste0("^", field, ":"), pdata[[c]], ignore.case = TRUE)), TRUE)]
  if (!length(hit)) stop("field not found: ", field)
  tolower(sub(paste0("^", field, ":\\s*"), "", pdata[[hit[1]]], ignore.case = TRUE))
}

derive_condition <- function(pdata, co) {
  if (co$cond_field == "diagnosis_blank") {
    # ERCB: living donors carry an empty 'diagnosis:' value; tumor
    # nephrectomies (matched by `control`) are controls as well.
    v <- get_field(pdata, "diagnosis")
    return(factor(ifelse(trimws(v) == "" | grepl(co$control, v, ignore.case = TRUE),
                         "Control", "Disease"),
                  levels = c("Control", "Disease")))
  }
  v <- get_field(pdata, co$cond_field)
  factor(ifelse(grepl(co$control, v, ignore.case = TRUE), "Control", "Disease"),
         levels = c("Control", "Disease"))
}

# Specific diagnosis ("indication") for the plan's Plot 4 heatmap. Applied to the
# same field derive_condition() reads; the first matching pattern wins, so the
# mixed and negated labels come before the single-disease ones. Controls keep
# the label "Control" whatever their wording. An unmatched Disease label is an
# error, so a new cohort cannot silently fall into the wrong row.
INDICATIONS <- c(
  "focal segmental glomerul[a-z]* sclerosis/minimal|glomerulosclerosis and minimal" = "FSGS/MCD",
  "diabetic|dkd|_dn$"                  = "DKD",
  "focal|fsgs"                         = "FSGS",
  "minimal change"                     = "MCD",
  "membranous"                         = "MN",
  "anca"                               = "ANCA vasculitis",
  "lupus"                              = "Lupus nephritis",
  "non-igan"                           = "Other GN",
  "iga|immunoglobulin a nephropathy"   = "IgAN",
  "hypertensive"                       = "Hypertensive nephropathy",
  "thin membrane"                      = "Thin membrane disease",
  "glomerular disease, other"          = "Other GN",
  "tumor nephrectom"                   = "Tumor nephrectomy",
  "renal biopsy specimens"             = "CKD (unspecified)"
)
derive_indication <- function(pdata, co, condition) {
  v <- get_field(pdata, if (co$cond_field == "diagnosis_blank") "diagnosis" else co$cond_field)
  out <- rep(NA_character_, length(v))
  for (pat in names(INDICATIONS))
    out[is.na(out) & grepl(pat, v, ignore.case = TRUE)] <- INDICATIONS[[pat]]
  out[condition == "Control"] <- "Control"
  if (anyNA(out)) stop(co$id, ": no indication for label(s): ",
                       paste(unique(v[is.na(out)]), collapse = " | "))
  data.table(Label = trimws(v), Indication = out)
}

derive_compartment <- function(pdata, co) {
  if (co$compartment != "split") return(rep(co$compartment, nrow(pdata)))
  # GSE30122 interleaves both compartments in one matrix.
  s <- tolower(pdata$source_name_ch1)
  ifelse(grepl("glomerul", s), "Glomerular",
         ifelse(grepl("tubul", s), "Tubulointerstitial", NA_character_))
}

# ---------------------------------------------------------------------------
# 3. Probe -> SYMBOL for the 5 panel genes, then collapse multiple probes per
#    gene with WGCNA::collapseRows(method = "MaxMean"): the probe with the
#    highest mean expression represents the gene.
# ---------------------------------------------------------------------------
panel_probes <- function(idtype, probes) {
  if (idtype == "symbol") {   # RNA-Seq deposits: rows are gene symbols already
    hit <- intersect(cfg$genes, probes)
    return(setNames(hit, hit))
  }
  if (idtype == "entrez") {   # "<entrez>_at", or bare "<entrez>" (GSE35488, CDF v10)
    want <- c(paste0(ENTREZ, "_at"), ENTREZ)
    m <- setNames(rep(names(ENTREZ), 2), want)
    hit <- intersect(want, probes)
    return(setNames(unname(m[hit]), hit))
  }
  if (!requireNamespace(idtype, quietly = TRUE))
    stop("annotation package not installed: ", idtype)
  # The ChipDb object is exported under the full package name, .db suffix included.
  db  <- getExportedValue(idtype, idtype)
  sym <- AnnotationDbi::mapIds(db, keys = probes, column = "SYMBOL",
                               keytype = "PROBEID", multiVals = "first")
  sym <- sym[!is.na(sym) & sym %in% cfg$genes]
  setNames(unname(sym), names(sym))
}

extract_panel <- function(sm, co) {
  expr <- sm$expr
  map  <- panel_probes(co$idtype, rownames(expr))
  if (!length(map)) return(NULL)
  sub  <- expr[names(map), , drop = FALSE]
  if (isTRUE(co$log2)) {
    stopifnot(min(sub, na.rm = TRUE) > -1)
    sub <- log2(sub + 1)                      # Agilent linear -> log2
  }
  # MaxMean collapse; rows of the result are named by gene symbol.
  cr <- WGCNA::collapseRows(sub, rowGroup = unname(map[rownames(sub)]),
                            rowID = rownames(sub), method = "MaxMean")
  cr$datETcollapsed
}

# ---------------------------------------------------------------------------
# 4. Load every cohort
# ---------------------------------------------------------------------------
panels <- list(); phenos <- list(); report <- list(); indic <- list()
for (co in COHORTS) {
  f <- list.files(file.path(cfg$raw_dir, co$dir), "_series_matrix\\.txt$", full.names = TRUE)
  f <- f[grepl(co$id, basename(f), fixed = TRUE)]
  if (!length(f)) { warning(co$id, ": series matrix not found"); next }

  sm  <- read_sm(f[1])
  if (!is.null(co$expr_file)) {   # RNA-Seq: the series-matrix table is empty
    ex <- data.table::fread(file.path(cfg$raw_dir, co$dir, co$expr_file))
    sm$expr <- as.matrix(ex[, -1]); rownames(sm$expr) <- ex[[1]]
    # The file may omit samples (GSE175759: submitter-flagged technical outliers);
    # series-matrix samples without values are dropped.
    stopifnot(all(colnames(sm$expr) %in% rownames(sm$pheno)))
    sm$pheno <- sm$pheno[rownames(sm$pheno) %in% colnames(sm$expr), , drop = FALSE]
    sm$expr <- sm$expr[, rownames(sm$pheno)]
  }
  mat <- tryCatch(extract_panel(sm, co),
                  error = function(e) { warning(co$id, ": ", conditionMessage(e)); NULL })
  if (is.null(mat)) next

  ph <- data.table(GSM = colnames(mat), StudyID = co$id,
                   Condition   = derive_condition(sm$pheno, co),
                   Compartment = derive_compartment(sm$pheno, co))
  ph <- cbind(ph, derive_indication(sm$pheno, co, ph$Condition))
  ph <- ph[!is.na(Compartment)]
  mat <- mat[, ph$GSM, drop = FALSE]
  # Indication travels in its own table: the Stage 4 / pooling input keeps its
  # fixed 7-column layout.
  indic[[co$id]] <- ph[, .(StudyID, GSM, Compartment, Condition, Indication, Label)]
  ph[, c("Indication", "Label") := NULL]

  panels[[co$id]] <- mat; phenos[[co$id]] <- ph
  report[[co$id]] <- data.table(StudyID = co$id, genes = nrow(mat),
                                n = ncol(mat),
                                Control = sum(ph$Condition == "Control"),
                                Disease = sum(ph$Condition == "Disease"),
                                median  = round(median(mat, na.rm = TRUE), 2))
  msg(co$id, ": ", nrow(mat), " genes x ", ncol(mat), " samples")
}
stopifnot(length(panels) > 1)
qc <- rbindlist(report); print(qc)
write_xlsx(qc, file.path(cfg$out_dir, "cohort_qc.xlsx"))
indic <- rbindlist(indic)
write_xlsx(indic, file.path(cfg$out_dir, "sample_indication.xlsx"))
msg("indications (all cohorts):")
print(indic[, .N, by = .(Indication)][order(-N)])

# ---------------------------------------------------------------------------
# 5. Merge on the common gene set, then correct the aggregated samples (all
#    datasets except cfg$subsets) in ONE model.
#    removeBatchEffect fits each gene independently, so subsetting to the panel
#    before correction gives results identical to correcting the full matrix.
# ---------------------------------------------------------------------------
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
# GSE30122 shares all 22 GSM accessions with each of GSE30528/GSE30529, so column
# names are NOT unique. Index positionally — `expr[, pheno$GSM]` would silently
# resolve every duplicate to the first matching column.
stopifnot(identical(unname(colnames(expr)), pheno$GSM))
pheno[, col := .I]

# ONE model over the aggregated samples.
#   batch  = StudyID; design = intercept + Disease — plan §3, "~ StudyID +
#            Condition" (by instruction, 2026-11-04). The Batch column (study x
#            compartment unit) is still the unit of Stage 4 and the pooling.
# Compartment main effects are deliberately NOT in the design: compartment is
# nested in the batch units, so they are not estimable, and removeBatchEffect
# would silently zero the aliased batch coefficients (`beta[is.na(beta)] <- 0`)
# and leave those studies uncorrected. The rank check below guards that.
pheno[, Batch := if (uniqueN(Compartment) > 1) paste(StudyID, Compartment, sep = ":")
                 else StudyID, by = StudyID]
# Drop batch arms too small to estimate.
ok <- pheno[, .N, by = .(Batch, Condition)][N >= cfg$min_group_n]
ph <- pheno[paste(Batch, Condition) %in% ok[, paste(Batch, Condition)]]
if (nrow(ph) < nrow(pheno))
  msg("dropped ", nrow(pheno) - nrow(ph), " samples in arms with n < ", cfg$min_group_n)
stopifnot(uniqueN(ph$Batch) > 1, uniqueN(ph$Condition) == 2)

# All samples of every dataset BEFORE correction, with their batch unit — the
# input for Stage 4 and the pooling strategy, which model batch or compute
# effects within study. Single-dataset data, so the subsets are kept here.
unc <- melt(cbind(ph[, !"col"], t(expr[, ph$col, drop = FALSE])),
            id.vars = c("GSM", "StudyID", "Condition", "Compartment", "Batch"),
            variable.name = "Gene", value.name = "Expression")
write_xlsx(unc[!is.na(Expression)], file.path(cfg$out_dir, "panel_uncorrected_long.xlsx"))

# From here on everything is aggregated across datasets: leave out the subsets
# of GSE30122, whose samples it already contains.
ph <- ph[!StudyID %in% cfg$subsets]
stopifnot(!anyDuplicated(ph$GSM))
ph[, Compartment := factor(Compartment)]
batch <- factor(ph$StudyID)   # plan §3 Stage 3: ~ StudyID + Condition
# The same model, fitted gene by gene on the samples that have the gene: a batch
# unit with no values for a gene would otherwise make the design rank-deficient
# and removeBatchEffect would silently zero a batch coefficient. For a gene
# present everywhere this equals the all-genes call (it fits genes independently).
adj <- expr[, ph$col, drop = FALSE]
for (g in rownames(adj)) {
  ok <- !is.na(adj[g, ])
  b  <- droplevels(batch[ok])
  des <- model.matrix(~ Condition, data = ph[ok])
  X.batch <- model.matrix(~ b, contrasts.arg = list(b = contr.sum(levels(b))))[, -1, drop = FALSE]
  stopifnot("design + batch is rank-deficient: some batch effects would be silently zeroed" =
              qr(cbind(des, X.batch))$rank == ncol(des) + ncol(X.batch))
  adj[g, ok] <- limma::removeBatchEffect(adj[g, ok, drop = FALSE], batch = b, design = des)
}
msg("one model: corrected ", ncol(adj), " samples across ", nlevels(batch),
    " batch units = studies (", uniqueN(ph$StudyID), " studies, ", nlevels(ph$Compartment),
    " compartments)", if (length(cfg$subsets)) paste0("; not aggregated: ", paste(cfg$subsets, collapse = ", ")))
ph[, Compartment := as.character(Compartment)]

wide <- cbind(ph[, !"Batch"], t(adj))
wide[, col := NULL]        # positional index — not a measure variable
long <- melt(wide,
             id.vars = c("GSM", "StudyID", "Condition", "Compartment"),
             variable.name = "Gene", value.name = "Expression")
long <- long[!is.na(Expression)]
long[, Gene := factor(Gene, levels = genes)]
write_xlsx(long, file.path(cfg$out_dir, "panel_harmonised_long.xlsx"))

# Residual scale check. removeBatchEffect fixes LOCATION, not spread, so the
# question is whether per-study spread still differs AFTER correction. This must
# be measured per gene — pooling genes would be dominated by between-gene
# variance and would look falsely uniform.
sd_tab <- long[, .(sd = sd(Expression, na.rm = TRUE)), by = .(Compartment, Gene, StudyID)]
sd_rng <- sd_tab[, .(min_sd = round(min(sd), 2), max_sd = round(max(sd), 2),
                     ratio  = round(max(sd) / min(sd), 1)),
                 by = .(Compartment, Gene)][order(-ratio)]
msg("per-gene spread across studies after correction (ratio >3 => also scale):")
print(sd_rng)
write_xlsx(sd_tab, file.path(cfg$out_dir, "residual_sd_by_study_gene.xlsx"))
write_xlsx(sd_rng, file.path(cfg$out_dir, "residual_sd_ratio.xlsx"))

# ---------------------------------------------------------------------------
# 6. Plot — 5 genes on one axis, Control vs Disease, faceted by compartment
# ---------------------------------------------------------------------------
pal <- c(Control = "#4DAF4A", Disease = "#E41A1C")

sub_lab <- sprintf("%d cohorts, %d samples", uniqueN(long$StudyID), uniqueN(long$GSM))
n_lab <- long[Gene == levels(long$Gene)[1], .N, by = .(Compartment, Condition)]

# --- gridline control, driven by cfg$gridlines -------------------------------
# Vertical rules are a geom_vline layer (on group BOUNDARIES, x = 1.5, 2.5, ...),
# not theme gridlines — a discrete-axis x gridline would fall at x = 1, 2, 3,
# i.e. straight through the middle of each pair instead of between them.
grid_theme <- function() {
  if (cfg$gridlines %in% c("horizontal", "both"))
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
  else
    theme(panel.grid = element_blank())
}
sep_layer <- function(d) {
  if (!cfg$gridlines %in% c("vertical", "both")) return(NULL)
  geom_vline(data = d, aes(xintercept = xintercept), inherit.aes = FALSE,
             colour = "grey78", linewidth = 0.3)
}

# Interior boundaries only; the compartment facets all share the same 5 genes.
sep_gene <- data.table(xintercept = seq_len(nlevels(long$Gene) - 1L) + 0.5)
# Each compartment facet holds a different number of cohorts (8 / 2 / 2), so
# these must be computed per facet rather than shared.
sep_study <- long[, .(n = uniqueN(StudyID)), by = Compartment][
  n > 1, .(xintercept = seq_len(n - 1L) + 0.5), by = Compartment]

p <- ggplot(long, aes(Gene, Expression, fill = Condition)) +
  sep_layer(sep_gene) +
  geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.45,
               position = position_dodge(0.78), alpha = 0.9) +
  # Points carry the per-sample distribution, but at n>600 they swamp the boxes;
  # thin the dense strata so the summary stays readable.
  geom_point(data = long[, .SD[if (.N > 120) sample(.N, 120) else seq_len(.N)],
                         by = .(Compartment, Gene, Condition)],
             position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.78),
             size = 0.55, alpha = 0.22, colour = "grey20") +
  facet_wrap(~ Compartment, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = pal) +
  labs(
    title    = "RenoRisk panel across curated GEO kidney-biopsy cohorts",
    subtitle = sub_lab,
    # Short label: at cfg$height = 3.25 the panel is ~1.5in tall and a long
    # rotated axis title clips. The qualifier lives in the subtitle and caption.
    x = NULL, y = "Expression (log2)",
    caption = "Submitter-processed series-matrix values; study offsets removed, scale differences not."
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 9),
        plot.caption = element_text(colour = "grey45", size = 7.5),
        axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
        axis.title.y = element_text(size = 10),
        legend.position = "top") +
  grid_theme()

out <- file.path(cfg$fig_dir, "panel_boxplot_by_compartment.png")
ggsave(out, p, width = cfg$width, height = cfg$height, dpi = cfg$dpi)
msg("wrote ", out)

# Per-cohort view — shows whether the correction actually aligned the studies.
p2 <- ggplot(long, aes(StudyID, Expression, fill = Condition)) +
  sep_layer(sep_study) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4, linewidth = 0.4,
               width = 0.36, position = position_dodge(0.75)) +
  facet_grid(Gene ~ Compartment, scales = "free", space = "free_x") +
  scale_fill_manual(values = pal) +
  labs(title = "Per-cohort check after batch correction",
       subtitle = sub_lab,
       x = NULL, y = "Batch-corrected expression (log2)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        # extra left margin: the rotated first cohort label overruns the panel
        plot.margin = margin(5.5, 5.5, 5.5, 34),
        legend.position = "top") +
  grid_theme()
ggsave(file.path(cfg$fig_dir, "panel_boxplot_per_cohort.png"), p2,
       width = cfg$pc_width, height = cfg$pc_height, dpi = cfg$dpi)

# ---------------------------------------------------------------------------
# 7. PCA before vs after batch correction (plan §7, Plot 1)
#    Same samples, same 5 genes, PCA fitted separately on each matrix. Genes are
#    scaled to unit variance: the per-gene-centred cohorts make raw variances
#    incomparable, and unscaled PCA would be driven by the highest-variance gene.
#    The R² table quantifies what the plot shows: the share of each PC explained
#    by batch unit (should fall after correction) and by Condition (should stay).
# ---------------------------------------------------------------------------
pca_one <- function(m, when) {
  pc  <- prcomp(t(m), center = TRUE, scale. = TRUE)
  ve  <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  sc  <- data.table(ph[, .(StudyID, Batch, Compartment, Condition)],
                    PC1 = pc$x[, 1], PC2 = pc$x[, 2], When = when)
  r2  <- rbindlist(lapply(c("PC1", "PC2"), function(k) data.table(
    When = when, PC = k, var_explained_pct = round(ve[as.integer(sub("PC", "", k))], 1),
    R2_batch     = round(summary(lm(sc[[k]] ~ sc$Batch))$r.squared, 3),
    R2_condition = round(summary(lm(sc[[k]] ~ sc$Condition))$r.squared, 3))))
  list(scores = sc, r2 = r2)
}
cg  <- intersect(cfg$genes, common)   # PCA needs complete rows: genes on every cohort
pcs <- list(pca_one(expr[cg, ph$col, drop = FALSE], "Before correction"),
            pca_one(adj[cg, , drop = FALSE], "After correction"))
pca_r2 <- rbindlist(lapply(pcs, `[[`, "r2"))
write_xlsx(pca_r2, file.path(cfg$out_dir, "pca_batch_r2.xlsx"))
msg("PCA, share of PC variance explained by batch unit / Condition:")
print(pca_r2)

scores <- rbindlist(lapply(pcs, `[[`, "scores"))
strip <- pca_r2[, .(lab = sprintf("%s\nPC1 %.0f%%, PC2 %.0f%% of variance; batch R² PC1 %.2f, PC2 %.2f",
                                  When[1], var_explained_pct[1], var_explained_pct[2],
                                  R2_batch[1], R2_batch[2])), by = When]
scores[strip, When := i.lab, on = "When"]
scores[, When := factor(When, levels = strip$lab)]

p3 <- ggplot(scores, aes(PC1, PC2, colour = StudyID, shape = Compartment)) +
  geom_point(size = 1.3, alpha = 0.75) +
  facet_wrap(~ When, nrow = 1, scales = "free") +
  scale_shape_manual(values = c(16, 17, 15)) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1), ncol = 1),
         shape  = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(title = sprintf("PCA of the panel (%s) before and after batch correction",
                       paste(cg, collapse = ", ")),
       subtitle = sub_lab,
       caption = paste("Genes centred and scaled to unit variance; PCA fitted separately per panel.",
                       sprintf("Batch R² = share of the PC explained by the %d study × compartment units.",
                               nlevels(batch)))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 9),
        plot.caption = element_text(colour = "grey45", size = 7.5),
        strip.text = element_text(size = 8.5))
ggsave(file.path(cfg$fig_dir, "pca_before_after_correction.png"), p3,
       width = cfg$pca_width, height = cfg$pca_height, dpi = cfg$dpi)
msg("wrote ", file.path(cfg$fig_dir, "pca_before_after_correction.png"))

# ---------------------------------------------------------------------------
# 8. UMAP before vs after batch correction (plan §7, Plot 1 — the UMAP half)
#    Same scaled 5-gene matrices as the PCA. UMAP distances are not linear, so
#    an R² does not apply; batch mixing is summarised instead by the mean
#    silhouette width of the batch units on the UMAP coordinates (near 1 =
#    studies form separate clusters, near 0 or below = studies are mixed).
#    Single-threaded with a fixed seed so the layout is reproducible.
# ---------------------------------------------------------------------------
umap_one <- function(m, when) {
  xy <- uwot::umap(scale(t(m)), n_neighbors = cfg$umap_neighbors,
                   min_dist = cfg$umap_min_dist, n_threads = 1, seed = 1)
  sil <- mean(cluster::silhouette(as.integer(factor(ph$Batch)), dist(xy))[, "sil_width"])
  list(coords = data.table(ph[, .(StudyID, Batch, Compartment, Condition)],
                           UMAP1 = xy[, 1], UMAP2 = xy[, 2], When = when),
       sil = data.table(When = when, batch_silhouette = round(sil, 3)))
}
ums <- list(umap_one(expr[cg, ph$col, drop = FALSE], "Before correction"),
            umap_one(adj[cg, , drop = FALSE], "After correction"))
umap_sil <- rbindlist(lapply(ums, `[[`, "sil"))
write_xlsx(umap_sil, file.path(cfg$out_dir, "umap_batch_silhouette.xlsx"))
msg("UMAP, mean silhouette width of the batch units:")
print(umap_sil)

ucoords <- rbindlist(lapply(ums, `[[`, "coords"))
ustrip <- umap_sil[, .(lab = sprintf("%s\nbatch silhouette %.2f", When, batch_silhouette)), by = When]
ucoords[ustrip, When := i.lab, on = "When"]
ucoords[, When := factor(When, levels = ustrip$lab)]

p4 <- ggplot(ucoords, aes(UMAP1, UMAP2, colour = StudyID, shape = Compartment)) +
  geom_point(size = 1.3, alpha = 0.75) +
  facet_wrap(~ When, nrow = 1, scales = "free") +
  scale_shape_manual(values = c(16, 17, 15)) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1), ncol = 1),
         shape  = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(title = sprintf("UMAP of the panel (%s) before and after batch correction",
                       paste(cg, collapse = ", ")),
       subtitle = sub_lab,
       caption = sprintf(paste("Genes centred and scaled to unit variance; UMAP fitted separately per panel",
                               "(n_neighbors = %d, min_dist = %.2f, seed 1). Batch silhouette = mean",
                               "silhouette width of the %d study × compartment units."),
                         cfg$umap_neighbors, cfg$umap_min_dist, nlevels(batch))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 9),
        plot.caption = element_text(colour = "grey45", size = 7.5),
        strip.text = element_text(size = 8.5))
ggsave(file.path(cfg$fig_dir, "umap_before_after_correction.png"), p4,
       width = cfg$pca_width, height = cfg$pca_height, dpi = cfg$dpi)
msg("wrote ", file.path(cfg$fig_dir, "umap_before_after_correction.png"))

writeLines(capture.output(sessionInfo()),
           file.path(cfg$out_dir, "stage23_sessionInfo.txt"))
msg("done")
