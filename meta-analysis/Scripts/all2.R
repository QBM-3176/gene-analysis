#!/usr/bin/env Rscript
# =============================================================================
# RenoRisk: sample-level ML input with every gene common to all series.
#
# Self-contained: reads only Data/Raw_filtered inputs, no other script's output.
# Builds the rows all.R Part I builds for all.xlsx (931 microarray samples, same metadata,
# Sex, Disease_Label) with the same rules - condition / compartment /
# indication per the cohort registry, batch unit = study x compartment, arms
# with n < 3 dropped, sex reported else inferred from XIST / Y genes - and
# carries every gene measured in all 19 series (14 microarray + 5 RNA-Seq) instead of
# the five panel genes. Expression as all.R builds the panel: series-matrix values as
# deposited (before cross-study batch correction), probes mapped to gene
# symbols, GSE66494 log2(x + 1), multi-probe genes collapsed with
# WGCNA::collapseRows(method = "MaxMean"). GSE30528/GSE30529 are not read:
# their samples are all in GSE30122.
#
# RNA-Seq (user decision, 2026-11-03): GSE142025, GSE141295, GSE175759,
# GSE137570-GPL23227, GSE219185 join as rows with Platform = "RNA-Seq". Their series
# matrices are metadata only; expression comes from <id>_expression.txt (GSM
# columns): submitter-normalised log2 (GSE219185: voom log2-CPM, Ensembl IDs ->
# symbol via org.Hs.eg.db), or DESeq2 VST made from the deposited counts by
# rnaseq_counts_vst.R (GSE141295, GSE175759), on a different scale from the arrays;
# no cross-platform correction is applied. GSE137570 and GSE219185 have no
# controls, so study identity predicts the label: only grouped-by-study CV
# (LODO) is valid. GSE219185 (NEPTUNE, tubulointerstitium) shares patients with
# GSE108109 (NEPTUNE, glomeruli): put both in one LODO group.
#
# Input : Data/Raw_filtered/<GSE>/<id>_series_matrix.txt (19 series matrices)
#         + <id>_expression.txt for the 5 RNA-Seq series
# Output: Data/ML_input/all2.xlsx; Results/ML_training/all2_genes.xlsx, all2_sex_by_cohort.xlsx,
#         all2_sessionInfo.txt
# Usage : Rscript Scripts/all2.R   (no arguments, no network; ~2 min)
# =============================================================================
set.seed(1)
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(writexl)
})

# 1. Configuration.
cfg <- list(
  project_dir = path.expand("~/Striki/Quant/quant.2/MetaAnalysis"),
  genes       = c("C4BPA", "COL13A1", "HABP2", "CFH", "CFI"),   # panel; reported if not common
  min_group_n = 3,           # drop a batch arm with fewer samples than this (all.R Stage 3)
  x_gene      = "XIST",
  y_genes     = c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY"),
  min_sex_sep = 4            # inferred sex used only if the two clusters are this
                             # many within-cluster SDs apart; else left unknown
)
cfg$raw_dir <- file.path(cfg$project_dir, "Data", "Raw_filtered")
cfg$out_dir <- file.path(cfg$project_dir, "Results", "ML_training")
cfg$ml_input_dir <- file.path(cfg$project_dir, "Data", "ML_input")   # all2.xlsx: ML input
for (d in c(cfg$out_dir, cfg$ml_input_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# The cohort registry: the all.R fields this script needs, for the 11
# aggregated series (all.R's SUBSETS GSE30528/GSE30529 left out).
#   idtype        probe-ID route: "entrez" (custom CDF) or a ChipDb package
#   compartment   tissue, or "split" (GSE30122 holds both, read per sample)
#   cond_field    where the condition label lives; `control` is a regex -> Control
#   blank_control ERCB: an EMPTY diagnosis is a living donor (Control)
#   log2          deposited on a linear scale (Agilent) -> log2(x + 1)
#   expr_file     RNA-Seq: expression table beside the (empty) series matrix;
#                 idtype "symbol" (rows are symbols) or "ensembl"; no `control`
#                 = the series has no control samples
COHORTS <- list(
  list(id = "GSE104948-GPL22945", dir = "GSE104948", idtype = "entrez",
       compartment = "Glomerular", cond_field = "diagnosis", blank_control = TRUE,
       control = "tumor nephrectom"),
  list(id = "GSE104954-GPL22945", dir = "GSE104954", idtype = "entrez",
       compartment = "Tubulointerstitial", cond_field = "diagnosis", blank_control = TRUE,
       control = "tumor nephrectom"),
  list(id = "GSE108109", dir = "GSE108109", idtype = "entrez",
       compartment = "Glomerular", cond_field = "source", control = "living donor"),
  list(id = "GSE200828", dir = "GSE200828", idtype = "entrez",
       compartment = "Glomerular", cond_field = "source", control = "living donor"),
  list(id = "GSE99339-GPL19184", dir = "GSE99339", idtype = "entrez",
       compartment = "Glomerular", cond_field = "description", control = "tumor nephrectomy"),
  # Added 2026-11-03 (by instruction); GSE35488 has bare Entrez IDs (CDF v10).
  list(id = "GSE93798", dir = "GSE93798", idtype = "entrez",
       compartment = "Glomerular", cond_field = "source", control = "control"),
  list(id = "GSE35488", dir = "GSE35488", idtype = "entrez",
       compartment = "Tubulointerstitial", cond_field = "source", control = "living donor"),
  list(id = "GSE111154", dir = "GSE111154", idtype = "hta20transcriptcluster.db",
       compartment = "Whole biopsy", cond_field = "source", control = "control"),
  list(id = "GSE121233", dir = "GSE121233", idtype = "hta20transcriptcluster.db",
       compartment = "Glomerular", cond_field = "disease state", control = "^control"),
  list(id = "GSE129973", dir = "GSE129973", idtype = "hta20transcriptcluster.db",
       compartment = "Glomerular", cond_field = "source", control = "tumor nephrectom"),
  list(id = "GSE96804", dir = "GSE96804", idtype = "hta20transcriptcluster.db",
       compartment = "Glomerular", cond_field = "source", control = "tumor nephrectom"),
  list(id = "GSE30122", dir = "GSE30122", idtype = "hgu133a2.db",
       compartment = "split", cond_field = "source", control = "control kidney"),
  list(id = "GSE116626", dir = "GSE116626", idtype = "illuminaHumanv4.db",
       compartment = "Whole biopsy", cond_field = "subject status", control = "living donor"),
  list(id = "GSE66494", dir = "GSE66494", idtype = "hgug4112a.db", log2 = TRUE,
       compartment = "Whole biopsy", cond_field = "sample type", control = "normal kidney"),
  list(id = "GSE142025", dir = "GSE142025", idtype = "symbol", platform = "RNA-Seq",
       expr_file = "GSE142025_expression.txt",
       compartment = "Whole biopsy", cond_field = "group", control = "^control$"),
  # RNA-Seq from raw counts (DESeq2 VST by rnaseq_counts_vst.R).
  list(id = "GSE141295", dir = "GSE141295", idtype = "symbol", platform = "RNA-Seq",
       expr_file = "GSE141295_expression.txt",
       compartment = "Glomerular", cond_field = "tissue1", control = "normal cortex"),
  list(id = "GSE175759", dir = "GSE175759", idtype = "symbol", platform = "RNA-Seq",
       expr_file = "GSE175759_expression.txt",
       compartment = "Tubulointerstitial", cond_field = "diagnosis", control = "^control$"),
  list(id = "GSE137570-GPL23227", dir = "GSE137570", idtype = "symbol", platform = "RNA-Seq",
       expr_file = "GSE137570-GPL23227_expression.txt",
       compartment = "Whole biopsy", cond_field = "tissue"),
  list(id = "GSE219185", dir = "GSE219185", idtype = "ensembl", platform = "RNA-Seq",
       expr_file = "GSE219185_expression.txt",
       compartment = "Tubulointerstitial", cond_field = "disease")
)
names(COHORTS) <- vapply(COHORTS, `[[`, "", "id")

# Specific diagnosis per sample, matched on the condition field; first match
# wins. Controls are "Control"; an unmatched Disease label is an error.
INDICATIONS <- c(
  "focal segmental glomerul[a-z]* sclerosis/minimal|glomerulosclerosis and minimal" = "FSGS/MCD",
  "diabetic|dkd|^(early|advanced)_dn$" = "DKD",
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
  "renal biopsy specimens"    = "CKD (unspecified)",
  "^kidney$"                  = "CKD (unspecified)"    # GSE137570: tissue only, all CKD
)

# 2. Helpers (as in all.R Part B).
read_series_matrix <- function(path, expr_file = NULL) {
  lines <- readLines(path, warn = FALSE)
  meta  <- grep("^!Sample_", lines, value = TRUE)
  parts <- lapply(meta, function(ln) {
    p <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    list(key = sub("^!Sample_", "", p[1]), vals = gsub('^"|"$', "", p[-1]))
  })
  keys  <- make.unique(vapply(parts, `[[`, "", "key"))
  pheno <- as.data.frame(setNames(lapply(parts, `[[`, "vals"), keys), stringsAsFactors = FALSE)
  rownames(pheno) <- pheno$geo_accession
  if (!is.null(expr_file)) {                       # RNA-Seq: series matrix is metadata only
    dt <- fread(expr_file, sep = "\t", header = TRUE)
    expr <- as.matrix(dt[, -1L]); rownames(expr) <- dt[[1L]]
    # The file may omit samples (GSE175759: submitter-flagged technical outliers).
    stopifnot(all(colnames(expr) %in% pheno$geo_accession), !anyDuplicated(rownames(expr)))
    pheno <- pheno[pheno$geo_accession %in% colnames(expr), , drop = FALSE]
    return(list(pheno = pheno, expr = expr[, pheno$geo_accession, drop = FALSE]))
  }
  b <- which(trimws(lines) == "!series_matrix_table_begin") + 1
  e <- which(trimws(lines) == "!series_matrix_table_end")   - 1
  stopifnot(e >= b)
  expr <- as.matrix(read.table(text = paste(lines[b:e], collapse = "\n"), header = TRUE,
                               sep = "\t", row.names = 1, check.names = FALSE, quote = '"',
                               comment.char = ""))
  stopifnot(identical(colnames(expr), pheno$geo_accession))
  list(pheno = pheno, expr = expr)
}

# Raw values of a metadata field: "source", "description", or a
# characteristics key (the "key: " prefix removed).
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
  ctrl <- if (is.null(co$control)) rep(FALSE, length(v)) else grepl(co$control, v, ignore.case = TRUE)
  if (isTRUE(co$blank_control)) ctrl <- ctrl | trimws(v) == ""
  ifelse(ctrl, "Control", "Disease")
}

derive_indication <- function(pheno, co, condition) {
  v <- tolower(field_values(pheno, co$cond_field))
  out <- rep(NA_character_, length(v))
  for (pat in names(INDICATIONS))
    out[is.na(out) & grepl(pat, v, ignore.case = TRUE)] <- INDICATIONS[[pat]]
  out[condition == "Control"] <- "Control"
  if (anyNA(out)) stop(co$id, ": no indication for label(s): ",
                       paste(unique(v[is.na(out)]), collapse = " | "))
  out
}

derive_compartment <- function(pheno, co) {
  if (co$compartment != "split") return(rep(co$compartment, nrow(pheno)))
  s <- tolower(pheno$source_name_ch1)       # GSE30122 interleaves both
  ifelse(grepl("glomerul", s), "Glomerular",
         ifelse(grepl("tubul", s), "Tubulointerstitial", NA_character_))
}

# Probe -> SYMBOL for the whole array. Custom-CDF probes are "<entrez>_at";
# ChipDb objects are exported under the full package name, ".db" included.
probe_symbols <- function(idtype, probes) {
  if (idtype == "symbol") return(setNames(probes, probes))
  if (idtype == "entrez") {
    probes <- grep("^[0-9]+(_at)?$", probes, value = TRUE)   # bare Entrez: GSE35488
    db <- org.Hs.eg.db::org.Hs.eg.db; keys <- sub("_at$", "", probes); keytype <- "ENTREZID"
  } else if (idtype == "ensembl") {
    db <- org.Hs.eg.db::org.Hs.eg.db; keys <- sub("\\..*$", "", probes); keytype <- "ENSEMBL"
  } else {
    db <- getExportedValue(idtype, idtype); keys <- probes; keytype <- "PROBEID"
  }
  sym <- suppressMessages(AnnotationDbi::mapIds(db, keys = keys, column = "SYMBOL",
                                                keytype = keytype, multiVals = "first"))
  sym <- setNames(unname(sym), probes)
  sym[!is.na(sym)]
}

# Multi-probe genes -> the probe with the highest mean (MaxMean); `map` is a
# named vector probe -> symbol. Rows of the result are named by symbol.
collapse_genes <- function(x, map)
  WGCNA::collapseRows(x[names(map), , drop = FALSE], rowGroup = unname(map),
                      rowID = names(map), method = "MaxMean")$datETcollapsed

# 3. Load each series: metadata per sample, expression on the all.R scale.
CK <- list()
for (co in COHORTS) {
  sm <- read_series_matrix(file.path(cfg$raw_dir, co$dir, paste0(co$id, "_series_matrix.txt")),
                           if (!is.null(co$expr_file)) file.path(cfg$raw_dir, co$dir, co$expr_file))
  if (isTRUE(co$log2)) {                           # Agilent linear -> log2
    stopifnot(min(sm$expr, na.rm = TRUE) > -1)
    sm$expr <- log2(sm$expr + 1)
  }
  condition <- derive_condition(sm$pheno, co)
  ph <- data.table(GSM = sm$pheno$geo_accession, StudyID = co$id,
                   Platform = if (is.null(co$platform)) "Microarray" else co$platform,
                   Compartment = derive_compartment(sm$pheno, co), Condition = condition,
                   Indication = derive_indication(sm$pheno, co, condition))
  CK[[co$id]] <- list(co = co, sm = sm, ph = ph[!is.na(Compartment)],
                      sym = probe_symbols(co$idtype, rownames(sm$expr)))
  message(co$id, ": ", nrow(sm$expr), " probes, ", length(unique(CK[[co$id]]$sym)),
          " genes, ", nrow(sm$pheno), " samples")
}

# 4. Samples: batch unit = study x compartment; arms below min_group_n dropped.
samples <- rbindlist(lapply(CK, `[[`, "ph"))
samples[, Batch := if (uniqueN(Compartment) > 1) paste(StudyID, Compartment, sep = ":")
                   else StudyID, by = StudyID]
samples[, n_arm := .N, by = .(Batch, Condition)]
if (any(samples$n_arm < cfg$min_group_n))
  message("dropped ", sum(samples$n_arm < cfg$min_group_n), " samples in arms with n < ",
          cfg$min_group_n)
samples <- samples[n_arm >= cfg$min_group_n][, n_arm := NULL]

# 5. Sex: reported, else inferred. Score = mean z(Y genes) - z(XIST) within
#    the study's samples, split by 2-means; accepted only if well separated.
reported_sex <- function(pheno) {
  v <- tryCatch(field_values(pheno, "sex"), error = function(e) NULL)
  if (is.null(v)) return(rep(NA_character_, nrow(pheno)))
  v <- tolower(trimws(v))
  ifelse(v %in% c("m", "male"), "M", ifelse(v %in% c("f", "female"), "F", NA_character_))
}

infer_sex <- function(ck, gsm) {
  map <- ck$sym[ck$sym %in% c(cfg$x_gene, cfg$y_genes)]
  y <- intersect(cfg$y_genes, map); has_x <- cfg$x_gene %in% map
  if (!has_x && !length(y))
    return(list(sex = rep(NA_character_, length(gsm)), sep = NA_real_,
                note = "sex marker genes not on platform"))
  ex <- ck$sm$expr
  z  <- t(scale(t(collapse_genes(ex[, colnames(ex) %in% gsm, drop = FALSE], map))))
  score <- (if (length(y)) colMeans(z[y, , drop = FALSE]) else 0) -
           (if (has_x) z[cfg$x_gene, ] else 0)
  km  <- kmeans(score, centers = range(score))
  sep <- abs(diff(km$centers[, 1])) /
         sqrt(mean(tapply(score, km$cluster, function(s) if (length(s) > 1) var(s) else 0)))
  if (!is.finite(sep) || sep < cfg$min_sex_sep)
    return(list(sex = setNames(rep(NA_character_, length(score)), names(score)), sep = sep,
                note = "no clear XIST/Y bimodality"))
  list(sex = setNames(ifelse(km$cluster == which.max(km$centers[, 1]), "M", "F"), names(score)),
       sep = sep, note = paste(c(if (length(y)) sprintf("%d Y genes", length(y)),
                                 if (has_x) "XIST"), collapse = " + "))
}

sex_tab <- list(); sex_log <- list()
for (sid in names(CK)) {
  ck  <- CK[[sid]]
  gsm <- colnames(ck$sm$expr)[colnames(ck$sm$expr) %in% samples[StudyID == sid, GSM]]
  rep_sex <- reported_sex(ck$sm$pheno)
  if (sid == "GSE93798")   # sex is in the title only: "... H12_M_67_1_Con_L"
    rep_sex <- sub("^.* [A-Z]+[0-9]+_([MF])_.*$", "\\1", ck$sm$pheno$title)
  rep_sex <- setNames(rep_sex, ck$sm$pheno$geo_accession)[gsm]
  inf <- infer_sex(ck, gsm)
  sex_tab[[sid]] <- data.table(StudyID = sid, GSM = gsm, Reported = unname(rep_sex),
                               Inferred = unname(inf$sex[gsm]))
  sex_log[[sid]] <- data.table(StudyID = sid, n = length(gsm), separation = round(inf$sep, 1),
                               inference = inf$note, reported = sum(!is.na(rep_sex)))
}
sex_tab <- rbindlist(sex_tab)
sex_tab[, Sex := fifelse(!is.na(Reported), Reported, Inferred)]
sex_log <- merge(rbindlist(sex_log),
                 sex_tab[, .(male = sum(Sex == "M", na.rm = TRUE),
                             female = sum(Sex == "F", na.rm = TRUE),
                             unknown = sum(is.na(Sex)),
                             concordant = sum(Reported == Inferred, na.rm = TRUE),
                             both_known = sum(!is.na(Reported) & !is.na(Inferred))),
                         by = StudyID], by = "StudyID")
samples <- merge(samples, sex_tab[, .(StudyID, GSM, Sex)], by = c("StudyID", "GSM"),
                 all.x = TRUE, sort = FALSE)
samples[, Disease_Label := as.integer(Condition == "Disease")]

# 6. Genes on every platform, collapsed MaxMean per series (per gene, so
#    restricting to the common genes first does not change the result).
meta_cols <- c("GSM", "StudyID", "Batch", "Compartment", "Condition",
               "Disease_Label", "Indication", "Sex", "Platform")
common <- sort(Reduce(intersect, lapply(CK, function(ck) unique(ck$sym))))
stopifnot(!any(common %in% meta_cols))
# C4BPA is absent from the GSE142025 and GSE137570 tables (not deposited), so it
# drops out of the intersection; reported, not fatal.
panel_missing <- setdiff(cfg$genes, common)
if (length(panel_missing)) message("panel genes not common to all series: ",
                                   paste(panel_missing, collapse = ", "))
message("genes common to all ", length(CK), " series: ", length(common))

expr_wide <- rbindlist(lapply(names(CK), function(id) {
  ck <- CK[[id]]
  gsm <- samples[StudyID == id, GSM]
  map <- ck$sym[ck$sym %in% common]
  x <- collapse_genes(ck$sm$expr[, gsm, drop = FALSE], map)
  x <- t(x[common, , drop = FALSE])
  stopifnot(all(is.finite(x)))
  data.table(StudyID = id, GSM = rownames(x), x)
}), use.names = TRUE)
setnames(expr_wide, c("StudyID", "GSM", common))

training_data <- merge(samples[, ..meta_cols], expr_wide, by = c("StudyID", "GSM"), sort = FALSE)
setcolorder(training_data, c(meta_cols, common))
setorder(training_data, StudyID, GSM)
stopifnot(
  nrow(training_data) == 931L + 36L + 24L + 87L + 24L + 220L,
  !anyDuplicated(training_data$GSM),
  uniqueN(training_data$StudyID) == 19L,
  uniqueN(training_data$Batch) == 20L,
  sum(training_data$Disease_Label == 1L) == 727L + 27L + 14L + 65L + 24L + 220L,
  sum(training_data$Disease_Label == 0L) == 204L + 9L + 10L + 22L,
  all(is.na(training_data$Sex) | training_data$Sex %in% c("M", "F")),
  all(vapply(training_data[, ..common], function(x) all(is.finite(x)), TRUE))
)

# 7. Gene coverage per series, and the README.
genes_tab <- data.table(
  Series = names(CK),
  Platform = vapply(CK, function(ck) ck$co$idtype, ""),
  Probes = vapply(CK, function(ck) nrow(ck$sm$expr), 0L),
  Genes_on_platform = vapply(CK, function(ck) length(unique(ck$sym)), 0L),
  Samples = vapply(names(CK), function(id) sum(training_data$StudyID == id), 0L),
  Common_genes = length(common))

readme <- data.table(
  Item = c("Rows", "Predictors", "Gene mapping", "Sex", "Target", "Metadata",
           "Expression scale", "Batch correction", "Excluded series",
           "Cross-study overlap", "Validation", "Input provenance"),
  Description = c(
    sprintf("%d unique GEO GSM accessions: %d disease and %d control records (the 931 all.xlsx microarray rows + 391 RNA-Seq rows: GSE142025 36, GSE141295 24, GSE175759 87, GSE137570-GPL23227 24, GSE219185 220; Platform column). Batch arms with n < %d dropped; batch unit = study x compartment (GSE30122 split).",
            nrow(training_data), sum(training_data$Disease_Label == 1L),
            sum(training_data$Disease_Label == 0L), cfg$min_group_n),
    sprintf("%d genes measured in all 19 series (columns after Platform, alphabetical), including %s%s; Sex is an optional categorical predictor.",
            length(common), paste(intersect(cfg$genes, common), collapse = ", "),
            if (length(panel_missing)) paste0(". Panel gene(s) NOT included (absent from at least one RNA-Seq table): ",
                                              paste(panel_missing, collapse = ", ")) else ""),
    "Probe -> gene symbol: org.Hs.eg.db for the Entrez custom-CDF cohorts (<entrez>_at; bare <entrez> for GSE35488), the platform ChipDb otherwise (first symbol per probe). Multi-probe genes collapsed with WGCNA::collapseRows (MaxMean), as in all.R. RNA-Seq: GSE142025/GSE137570 rows are symbols as deposited; GSE219185 Ensembl IDs -> symbol via org.Hs.eg.db, duplicates collapsed MaxMean. Unannotated probes are dropped. Matching is by symbol, so a gene with different symbols across annotation packages can be missed. Per-series coverage: all2_genes.xlsx.",
    "Sex: M = Male; F = Female; missing values remain NA. Reported sex takes priority; otherwise inferred from XIST vs Y-chromosome genes (2-means, accepted if clusters are >= 4 within-cluster SDs apart). Per-study summary: all2_sex_by_cohort.xlsx. Recorded/inferred sex, not gender identity.",
    "Disease_Label: 0 = Control; 1 = Disease. Condition preserves the source label.",
    "GSM, StudyID, Batch, Compartment and Indication describe samples; they are not model predictors.",
    "Series-matrix values as deposited (log2/relative expression); GSE66494 uses log2(x + 1). GSE30122 and GSE116626 are per-gene centred: they retain a relative scale for every gene. RNA-Seq rows (Platform = RNA-Seq) are normalised log2 from different pipelines (GSE219185 voom log2-CPM; GSE141295/GSE175759 DESeq2 VST computed from deposited counts), on a different scale from the arrays.",
    "No cross-study or cross-platform batch correction or additional ML scaling is applied.",
    "GSE30528 and GSE30529 are not read: all their samples are in GSE30122. GSE35487 (same arrays as GSE35488) and GSE37460 (arrays already in GSE104948/GSE99339) are not used.",
    "The expression audit flags GSE108109/GSE200828 as strongly overlapping and GSE121233/GSE129973 as additional candidates.",
    "Only grouped-by-study CV (LODO) is valid: GSE137570-GPL23227 and GSE219185 have no controls, so study identity predicts the label. GSE219185 (NEPTUNE tubulointerstitium) shares patients with GSE108109 (NEPTUNE glomeruli): put them in one LODO group. Use study/patient-aware validation. Keep linked studies together as appropriate; fit learned preprocessing and any gene selection only within training folds. No fixed train/test split is assigned here.",
    "Data/Raw_filtered series matrices (submitter-processed microarrays) and the five RNA-Seq <id>_expression.txt tables; no other pipeline output is read."
  ))

# 8. Write, then read back: the workbook must preserve values and labels.
output <- file.path(cfg$ml_input_dir, "all2.xlsx")
write_xlsx(list(training_data = training_data, README = readme), output)
write_xlsx(genes_tab, file.path(cfg$out_dir, "all2_genes.xlsx"))
write_xlsx(sex_log, file.path(cfg$out_dir, "all2_sex_by_cohort.xlsx"))
reloaded <- as.data.table(read_excel(output, sheet = "training_data",
  col_types = c(rep("text", 5L), "numeric", "text", "text", "text", rep("numeric", length(common)))))
# Compared to 1e-12, not identical(): an in-memory double can lose its last
# bit on the first xlsx round trip (relative error ~5e-16).
stopifnot(
  identical(names(reloaded), names(training_data)),
  identical(reloaded$GSM, training_data$GSM),
  isTRUE(all.equal(as.matrix(reloaded[, ..common]), as.matrix(training_data[, ..common]),
                   tolerance = 1e-12)),
  identical(as.integer(reloaded$Disease_Label), training_data$Disease_Label),
  identical(reloaded$Sex, training_data$Sex),
  identical(reloaded$Platform, training_data$Platform)
)
capture.output(sessionInfo(), file = file.path(cfg$out_dir, "all2_sessionInfo.txt"))
print(training_data[, .N, by = .(Condition)])
print(training_data[, .N, by = .(Sex)])
cat("Created", output, ":", nrow(training_data), "samples x", length(common), "genes\n")
