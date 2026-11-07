#!/usr/bin/env Rscript
# =============================================================================
# MetaAnalysis — RNA-Seq series WITHOUT controls: panel values for the ML inputs only
#
# By instruction (2026-11-03) RNA-Seq is in the same bag as the microarrays:
#   * GSE142025 (DN vs tumor-nephrectomy controls) is a regular Stage 2/3 cohort
#     (stage23_panel_boxplots.R) and so enters every analysis.
#   * GSE137570 (GPL23227; CKD patients only) and GSE219185 (NEPTUNE
#     tubulointerstitium; FSGS vs MCD) have NO control arm. Their batch offset
#     cannot be separated from disease, and they give no within-study
#     Disease-vs-Control effect, so they are kept out of batch correction, Stage 4,
#     pooling and the mega-analysis (user decision, 2026-11-03) and enter the ML
#     inputs only (create_ML_train_data.R, all.R Part I, all2.R).
#     The older GSE137570 GPL11154 part was removed from Raw_filtered by the user.
#
# Input : Data/Raw_filtered/<GSE>/ — the series matrix (labels; its table is
#         empty for RNA-Seq) and <id>_expression.txt (submitter-normalised log2,
#         gene x GSM; GSE219185 rows are Ensembl IDs, voom log2-CPM).
# Output: Results/RNASeq/
#           panel_rnaseq_long.xlsx  same 7 columns as Stage23 panel_uncorrected_long
#           rnaseq_samples.xlsx     StudyID, GSM, Compartment, Condition, Indication,
#                                   Label, Sex (reported, else inferred from XIST / Y genes)
#           rnaseq_qc.xlsx          one row per series
# Run BEFORE create_ML_train_data.R.
# =============================================================================

cfg <- list(
  project_dir = path.expand("~/Striki/Quant/quant.2/MetaAnalysis"),
  genes       = c("C4BPA", "COL13A1", "HABP2", "CFH", "CFI"),
  y_genes     = c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY"),
  x_gene      = "XIST",
  min_sex_sep = 4      # as Stage 4: accept inferred sex only if clusters are this far apart
)
cfg$raw_dir <- file.path(cfg$project_dir, "Data", "Raw_filtered")
cfg$out_dir <- file.path(cfg$project_dir, "Results", "RNASeq")

suppressPackageStartupMessages({ library(data.table); library(writexl) })
set.seed(1)
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

# cond_field = characteristic holding the diagnosis; indication maps its value
# (NULL: one indication for the whole series). No series here has controls.
COHORTS <- list(
  list(id = "GSE137570-GPL23227", gse = "GSE137570", idtype = "symbol",
       cond_field = "tissue", compartment = "Whole biopsy",
       indication = c(".*" = "CKD (unspecified)")),
  list(id = "GSE219185",          gse = "GSE219185", idtype = "ensembl",
       cond_field = "disease", compartment = "Tubulointerstitial",
       indication = c("focal segmental" = "FSGS", "minimal change" = "MCD")))

read_pheno <- function(path) {
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
infer_sex <- function(ex) {
  y <- intersect(cfg$y_genes, rownames(ex)); has_x <- cfg$x_gene %in% rownames(ex)
  if (!length(y) && !has_x) return(rep(NA_character_, ncol(ex)))
  z <- t(scale(t(ex[c(y, if (has_x) cfg$x_gene), , drop = FALSE])))
  score <- (if (length(y)) colMeans(z[y, , drop = FALSE]) else 0) -
           (if (has_x) z[cfg$x_gene, ] else 0)
  km  <- kmeans(score, centers = range(score))
  sep <- abs(diff(km$centers[, 1])) /
         sqrt(mean(tapply(score, km$cluster, function(s) if (length(s) > 1) var(s) else 0)))
  msg("  sex inference: ", length(y), " Y genes", if (has_x) " + XIST", ", separation ",
      round(sep, 1), " SD")
  if (!is.finite(sep) || sep < cfg$min_sex_sep) return(rep(NA_character_, ncol(ex)))
  ifelse(km$cluster == which.max(km$centers[, 1]), "M", "F")
}

long <- list(); smp <- list(); qc <- list()
for (co in COHORTS) {
  dir <- file.path(cfg$raw_dir, co$gse)
  ph  <- read_pheno(file.path(dir, paste0(co$id, "_series_matrix.txt")))
  lab <- ph$ch[[co$cond_field]]; stopifnot(!is.null(lab))
  ex  <- fread(file.path(dir, paste0(co$id, "_expression.txt")))
  ids <- ex[[1]]; ex <- as.matrix(ex[, -1])
  if (co$idtype == "ensembl") {   # Ensembl -> symbol, only for the genes used here
    sym <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = ids,
             column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
    ids <- unname(sym)
  }
  keep <- ids %in% c(cfg$genes, cfg$y_genes, cfg$x_gene)
  stopifnot(!anyDuplicated(ids[keep]))   # else a MaxMean collapse would be needed
  ex <- ex[keep, , drop = FALSE]; rownames(ex) <- ids[keep]
  stopifnot(setequal(colnames(ex), ph$gsm))
  ex <- ex[, ph$gsm, drop = FALSE]

  rep_sex <- ph$ch[["sex"]]
  rep_sex <- if (is.null(rep_sex)) rep(NA_character_, length(ph$gsm)) else
    c(male = "M", female = "F")[tolower(rep_sex)]
  sex <- ifelse(is.na(rep_sex), infer_sex(ex), rep_sex)
  ind <- rep(NA_character_, length(lab))
  for (pat in names(co$indication))
    ind[is.na(ind) & grepl(pat, lab, ignore.case = TRUE)] <- co$indication[[pat]]
  stopifnot(!anyNA(ind))

  smp[[co$id]] <- data.table(StudyID = co$id, GSM = ph$gsm, Compartment = co$compartment,
    Condition = "Disease", Indication = ind, Label = lab, Sex = unname(sex))
  genes <- intersect(cfg$genes, rownames(ex))
  long[[co$id]] <- data.table(GSM = rep(ph$gsm, each = length(genes)), StudyID = co$id,
    Condition = "Disease", Compartment = co$compartment, Batch = co$id,
    Gene = genes, Expression = as.numeric(ex[genes, ]))
  qc[[co$id]] <- data.table(Series = co$id, n = length(ph$gsm),
    Use = "ML inputs only (no control arm)",
    Indications = paste(names(table(ind)), table(ind), collapse = ", "),
    SexSource = if (all(!is.na(rep_sex))) "reported" else "inferred",
    PanelGenes = paste(genes, collapse = ", "),
    MissingPanelGenes = paste(setdiff(cfg$genes, genes), collapse = ", "))
}
long <- rbindlist(long); smp <- rbindlist(smp); qc <- rbindlist(qc)
stopifnot(nrow(long) > 0, !anyNA(long$Expression))

write_xlsx(long, file.path(cfg$out_dir, "panel_rnaseq_long.xlsx"))
write_xlsx(smp,  file.path(cfg$out_dir, "rnaseq_samples.xlsx"))
write_xlsx(qc,   file.path(cfg$out_dir, "rnaseq_qc.xlsx"))
print(qc)
print(smp[, .N, by = .(StudyID, Indication, Sex)])
writeLines(capture.output(sessionInfo()), file.path(cfg$out_dir, "rnaseq_sessionInfo.txt"))
msg("done: ", uniqueN(long$GSM), " samples")
