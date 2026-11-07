#!/usr/bin/env Rscript
# =============================================================================
# MetaAnalysis — RNA-Seq series deposited as RAW COUNTS only -> DESeq2 VST
#
# By instruction (2026-11-03): GSE175759 and GSE141295 are to be used like the
# other cohorts. Their series matrices are metadata-only and GEO holds the values
# only as per-sample count files in <id>_RAW.tar (downloaded into Data/Raw_filtered/
# <id>/ and extracted to <id>_RAW/).
# By instruction (2026-11-04): DESeq2 VST, the same transformation the per-dataset
# RNA-Seq boxplots use (Examples/boxplot_rr_dnaseq.r; plots_per_dataset.R
# process_rnaseq): design ~ condition, genes with rowSums >= 10, vst(blind = FALSE).
# Plan §3 Stage 3 (by instruction, 2026-11-04): "Apply ComBat-seq for raw integer
# RNA-Seq counts ... modeling ~ StudyID + Condition" — sva::ComBat_seq over both
# series on their common genes (batch = StudyID, group = Disease/Control) BEFORE
# the VST. Compartment is aliased with study here (glomeruli vs tubulointerstitium),
# so that difference is removed too; Stage 3 then corrects again with the arrays.
# Normalisation only — the submitters' alignment and counting are kept (no
# re-quantification; Stage 1 stays retired):
#   * Ensembl IDs (version stripped) -> symbol via org.Hs.eg.db; unmapped dropped;
#     counts of Ensembl IDs sharing a symbol are summed (counts are additive)
#   * condition = the cohort's label field (GSE141295 tissue1, GSE175759 diagnosis)
#   * GSE175759: the 3 samples the submitters flag "technical outlier:
#     Outlier_excluded from main analysis" are left out
#
#   GSE175759  featureCounts (header row), tubulointerstitium, 87 samples kept
#   GSE141295  HTSeq-count (no header; "__" summary rows dropped), glomeruli, 24
#
# Input : Data/Raw_filtered/<id>/<id>_RAW/GSM*_*.txt.gz, <id>_series_matrix.txt
# Output: Data/Raw_filtered/<id>/<id>_expression.txt (VST, gene symbol x GSM; read
#         by every script through the registry's `expr_file`);
#         Data/Raw_filtered/<id>/<id>_counts.txt (ComBat-seq counts; read by the
#         per-dataset DESeq2 boxplot branch);
#         Results/RNASeq/counts_vst_qc.xlsx, counts_vst_sessionInfo.txt
# Run BEFORE stage23_panel_boxplots.R (and all.R / all2.R).
# =============================================================================

cfg <- list(project_dir = path.expand("~/Striki/Quant/quant.2/MetaAnalysis"),
            genes = c("C4BPA", "COL13A1", "HABP2", "CFH", "CFI"),
            min_count = 10)
cfg$raw_dir <- file.path(cfg$project_dir, "Data", "Raw_filtered")
cfg$out_dir <- file.path(cfg$project_dir, "Results", "RNASeq")

suppressPackageStartupMessages({ library(data.table); library(DESeq2); library(writexl) })
set.seed(1)
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

SERIES <- list(
  list(id = "GSE175759", header = TRUE,  cond_field = "diagnosis", control = "^control$",
       drop_field = "technical outlier", drop = "outlier"),
  list(id = "GSE141295", header = FALSE, cond_field = "tissue1", control = "normal cortex"))

read_field <- function(lines, field) {   # characteristics value per sample
  ch <- grep("^!Sample_characteristics_ch1\t", lines, value = TRUE)
  v  <- lapply(strsplit(ch, "\t", fixed = TRUE), function(p) gsub('^"|"$', "", p[-1]))
  hit <- Filter(function(x) startsWith(tolower(x[1]), paste0(field, ":")), v)
  stopifnot(length(hit) == 1)
  trimws(sub("^[^:]*:", "", hit[[1]]))
}

raw <- list()
for (s in SERIES) {
  dir <- file.path(cfg$raw_dir, s$id)
  sm  <- readLines(file.path(dir, paste0(s$id, "_series_matrix.txt")), warn = FALSE)
  gsm <- gsub('"', "", strsplit(grep("^!Sample_geo_accession\t", sm, value = TRUE), "\t")[[1]][-1])
  files <- list.files(file.path(dir, paste0(s$id, "_RAW")), "^GSM[0-9]+_.*\\.txt\\.gz$", full.names = TRUE)
  names(files) <- sub("_.*$", "", basename(files))
  stopifnot(setequal(names(files), gsm))

  cnt <- lapply(files[gsm], function(f) {
    x <- fread(f, header = s$header, col.names = c("id", "count"))
    x[!startsWith(id, "__")]
  })
  ids <- cnt[[1]]$id
  stopifnot(all(vapply(cnt, function(x) identical(x$id, ids), TRUE)))
  counts <- do.call(cbind, lapply(cnt, `[[`, "count")); colnames(counts) <- gsm
  stopifnot(all(counts >= 0), all(counts == round(counts)))

  keep_smp <- rep(TRUE, length(gsm))
  if (!is.null(s$drop_field)) keep_smp <- !grepl(s$drop, read_field(sm, s$drop_field), ignore.case = TRUE)
  counts <- counts[, keep_smp, drop = FALSE]
  cond   <- factor(read_field(sm, s$cond_field)[keep_smp])

  sym <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
           keys = sub("\\..*$", "", ids), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
  ok  <- !is.na(sym)
  counts <- rowsum(counts[ok, , drop = FALSE], group = unname(sym[ok]), reorder = TRUE)
  raw[[s$id]] <- list(s = s, dir = dir, gsm = gsm, keep_smp = keep_smp, ids = ids,
                      counts = counts, cond = cond,
                      disease = ifelse(grepl(s$control, cond, ignore.case = TRUE), "Control", "Disease"))
}

# ComBat-seq on the genes both series measured.
common <- Reduce(intersect, lapply(raw, function(x) rownames(x$counts)))
stopifnot(all(cfg$genes %in% common))
all_cnt <- do.call(cbind, lapply(raw, function(x) x$counts[common, , drop = FALSE]))
study   <- rep(names(raw), vapply(raw, function(x) ncol(x$counts), 1L))
dis     <- unlist(lapply(raw, `[[`, "disease"), use.names = FALSE)
stopifnot(all(table(study, dis) > 0))
adj_cnt <- sva::ComBat_seq(all_cnt, batch = study, group = dis)
stopifnot(identical(dim(adj_cnt), dim(all_cnt)), all(adj_cnt >= 0))
msg("ComBat-seq: ", nrow(adj_cnt), " common genes x ", ncol(adj_cnt), " samples, ",
    length(raw), " studies")

qc <- list()
for (x in raw) {
  s <- x$s; dir <- x$dir; gsm <- x$gsm; keep_smp <- x$keep_smp; ids <- x$ids; cond <- x$cond
  counts <- adj_cnt[, study == s$id, drop = FALSE]
  fwrite(data.table(Gene = rownames(counts), counts),
         file.path(dir, paste0(s$id, "_counts.txt")), sep = "\t")

  # Same steps as the per-dataset boxplot branch, so both see identical values.
  dds <- DESeqDataSetFromMatrix(counts, data.frame(row.names = colnames(counts), condition = cond),
                                design = ~ condition)
  dds <- DESeq(dds[rowSums(counts(dds)) >= cfg$min_count, ], quiet = TRUE)
  v   <- SummarizedExperiment::assay(vst(dds, blind = FALSE))
  stopifnot(all(is.finite(v)))
  miss <- setdiff(cfg$genes, rownames(v))
  if (length(miss)) msg(s$id, ": below the count filter, absent from the VST table: ",
                        paste(miss, collapse = ", "))

  fwrite(data.table(Symbol = rownames(v), round(v, 6)),
         file.path(dir, paste0(s$id, "_expression.txt")), sep = "\t")
  pg <- intersect(cfg$genes, rownames(v))
  qc[[s$id]] <- data.table(Series = s$id, samples_deposited = length(gsm),
    samples_kept = ncol(counts), dropped = paste(gsm[!keep_smp], collapse = ", "),
    conditions = paste(sprintf("%s %d", levels(cond), tabulate(cond)), collapse = "; "),
    ensembl_ids = length(ids), symbols = nrow(x$counts), common_genes = nrow(counts), genes_vst = nrow(v),
    panel_missing = if (length(miss)) paste(miss, collapse = ", ") else NA_character_,
    size_factor_range = paste(signif(range(sizeFactors(dds)), 3), collapse = " - "),
    t(setNames(round(rowMeans(v[pg, , drop = FALSE]), 2), paste0("mean_", pg))))
  msg(s$id, ": ", ncol(v), " samples x ", nrow(v), " genes (VST) -> ",
      paste0(s$id, "_expression.txt"))
}
qc <- rbindlist(qc, fill = TRUE); print(qc)
write_xlsx(qc, file.path(cfg$out_dir, "counts_vst_qc.xlsx"))
writeLines(capture.output(sessionInfo()), file.path(cfg$out_dir, "counts_vst_sessionInfo.txt"))
