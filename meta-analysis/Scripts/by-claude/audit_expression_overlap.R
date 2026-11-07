#!/usr/bin/env Rscript
# Full-transcriptome overlap screen; does not alter the analysis or exclude data.
# Expression similarity is evidence of profile reuse, not proof of patient identity.
# Submitter processing, platforms, tissue and case mix limit sensitivity.
set.seed(1)
suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(writexl)
})
cfg <- list(project_dir = path.expand("~/Striki/Quant/quant.2/MetaAnalysis"),
            min_features = 1000L, screen_r = 0.95)
# Duplicate series (GSE30528, GSE30529, GSE35487, GSE37460; folders renamed -dup on
# 2026-11-03) are not read.
cfg$input_dir <- file.path(cfg$project_dir, "Data", "Raw_filtered")
cfg$out_dir <- file.path(cfg$project_dir, "Results", "SampleOverlap")
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
# Every series in use (2026-11-04: 19). `file` = basename prefix of the series
# matrix / expression file in Data/Raw_filtered/<study>/; RNA-Seq values come from
# <file>_expression.txt (symbol rows, or Ensembl for GSE219185).
registry <- data.table(
  study = c("GSE104948", "GSE104954", "GSE108109", "GSE200828", "GSE99339",
            "GSE93798", "GSE35488", "GSE121233", "GSE129973", "GSE96804", "GSE111154",
            "GSE30122", "GSE116626", "GSE66494",
            "GSE142025", "GSE141295", "GSE175759", "GSE137570", "GSE219185"),
  file = c("GSE104948-GPL22945", "GSE104954-GPL22945", "GSE108109", "GSE200828",
           "GSE99339-GPL19184", "GSE93798", "GSE35488", "GSE121233", "GSE129973",
           "GSE96804", "GSE111154", "GSE30122", "GSE116626", "GSE66494",
           "GSE142025", "GSE141295", "GSE175759", "GSE137570-GPL23227", "GSE219185"),
  namespace = c(rep("entrez", 7), rep("hta20transcriptcluster.db", 4),
                "hgu133a2.db", "illuminaHumanv4.db", "hgug4112a.db",
                rep("symbol", 4), "ensembl"))

# 1. Read the complete deposited matrix and sample metadata, without GEO fetches.
read_matrix <- function(study, file, namespace) {
  path <- file.path(cfg$input_dir, study, paste0(file, "_series_matrix.txt"))
  lines <- readLines(path, warn = FALSE)
  begin <- which(trimws(lines) == "!series_matrix_table_begin")
  end <- which(trimws(lines) == "!series_matrix_table_end")
  x <- if (namespace %in% c("symbol", "ensembl"))   # RNA-Seq: values in a separate file
    fread(file.path(cfg$input_dir, study, paste0(file, "_expression.txt")))
  else fread(text = paste(lines[(begin + 1L):(end - 1L)], collapse = "\n"),
             na.strings = c("NA", "null", ""))
  ids <- x[[1]]
  m <- as.matrix(x[, -1]); storage.mode(m) <- "double"
  rownames(m) <- ids
  metadata <- lapply(grep("^!Sample_", lines, value = TRUE), function(line) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    list(key = sub("^!Sample_", "", parts[1]),
         value = gsub('^"|"$', "", parts[-1]))
  })
  field <- function(key) {
    hit <- Filter(function(p) p$key == key, metadata)
    if (length(hit)) hit[[1]]$value else rep(NA_character_, ncol(m))
  }
  meta <- data.table(study = study, gsm = field("geo_accession"),
                     title = field("title"), source = field("source_name_ch1"))
  # An expression file may omit samples (GSE175759: submitter-flagged outliers).
  stopifnot(all(colnames(m) %in% meta$gsm))
  meta <- meta[match(colnames(m), gsm)]
  relations <- rbindlist(lapply(Filter(function(p) p$key == "relation", metadata),
    function(p) data.table(gsm = field("geo_accession"), relation = p$value)))
  if (nrow(relations)) relations <- relations[gsm %in% meta$gsm]
  if (study == "GSE66494") {
    m[m <= -1] <- NA_real_
    m <- log2(m + 1)
  }
  list(native = m, meta = meta, relations = relations)
}

# 2. Map to common Entrez IDs. Ambiguous mappings are excluded; mean collapse
#    is used for fingerprints only, not to replace the panel's MaxMean rule.
map_genes <- function(m, namespace) {
  if (namespace == "entrez") {   # "<entrez>_at", or bare "<entrez>" (GSE35488)
    ids <- ifelse(grepl("^[0-9]+(_at)?$", rownames(m)),
                  sub("_at$", "", rownames(m)), NA_character_)
  } else if (namespace %in% c("symbol", "ensembl")) {
    keys <- if (namespace == "ensembl") sub("\\..*$", "", rownames(m)) else rownames(m)
    ids <- unname(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys, column = "ENTREZID",
                                        keytype = toupper(namespace), multiVals = "asNA"))
  } else {
    db <- getExportedValue(namespace, namespace)
    keys <- intersect(rownames(m), AnnotationDbi::keys(db, keytype = "PROBEID"))
    mapped <- AnnotationDbi::mapIds(db, keys, column = "ENTREZID",
                                   keytype = "PROBEID", multiVals = "asNA")
    ids <- unname(mapped[rownames(m)])
  }
  keep <- !is.na(ids) & apply(m, 1L, function(v) all(is.finite(v)))
  ids <- ids[keep]; m <- m[keep, , drop = FALSE]
  total <- rowsum(m, group = ids, reorder = TRUE)
  total / as.numeric(table(factor(ids, levels = rownames(total))))
}
datasets <- list(); inventory <- list()
for (i in seq_len(nrow(registry))) {
  s <- registry$study[i]
  message("Reading and mapping ", s)
  d <- read_matrix(s, registry$file[i], registry$namespace[i])
  d$genes <- map_genes(d$native, registry$namespace[i])
  datasets[[s]] <- d
  inventory[[s]] <- data.table(study = s, samples = ncol(d$native),
    native_features = nrow(d$native), mapped_genes = nrow(d$genes),
    namespace = registry$namespace[i])
}
samples <- rbindlist(lapply(datasets, `[[`, "meta"))
relation_rows <- rbindlist(lapply(datasets, `[[`, "relations"), fill = TRUE)
relation_keys <- character()
if (nrow(relation_rows)) {
  for (i in seq_len(nrow(relation_rows))) {
    other <- regmatches(relation_rows$relation[i],
                        gregexpr("GSM[0-9]+", relation_rows$relation[i]))[[1]]
    for (o in other) relation_keys <- c(relation_keys,
      paste(sort(c(relation_rows$gsm[i], o)), collapse = "|"))
  }
}
relation_keys <- unique(relation_keys)

# 3. Profile correlations. Gene centering removes constitutive gene baselines.
#    Per-study gene standardization additionally reduces gene-wise scale shifts.
#    Neither uses disease labels; neither is training-data preprocessing.
normalize_columns <- function(m) {
  m <- sweep(m, 2L, colMeans(m), "-")
  norm <- sqrt(colSums(m^2)); norm[norm == 0] <- NA_real_
  sweep(m, 2L, norm, "/")
}
gene_center <- function(m) sweep(m, 1L, rowMeans(m), "-")
gene_scale <- function(m) {
  sd <- sqrt(rowSums(m^2) / (ncol(m) - 1L))
  # Drop numerically constant genes; no synthetic values or missing-value fill.
  sd[sd < 1e-8] <- NA_real_
  sweep(m, 1L, sd, "/")
}
similarity <- function(a, b) crossprod(normalize_columns(a), normalize_columns(b))
correlation_rows <- list(); summaries <- list(); k <- 0L
for (i in seq_along(datasets)) for (j in i:length(datasets)) {
  sa <- names(datasets)[i]; sb <- names(datasets)[j]
  da <- datasets[[i]]; db <- datasets[[j]]
  common <- sort(intersect(rownames(da$genes), rownames(db$genes)))
  a <- da$genes[common, , drop = FALSE]; b <- db$genes[common, , drop = FALSE]
  ca <- gene_center(a); cb <- gene_center(b)
  za <- gene_scale(ca); zb <- gene_scale(cb)
  good <- rowSums(!is.finite(za)) == 0L & rowSums(!is.finite(zb)) == 0L
  a <- a[good, , drop = FALSE]; b <- b[good, , drop = FALSE]
  ca <- ca[good, , drop = FALSE]; cb <- cb[good, , drop = FALSE]
  za <- za[good, , drop = FALSE]; zb <- zb[good, , drop = FALSE]
  if (nrow(a) < cfg$min_features) next
  message("Comparing ", sa, " / ", sb, " (", nrow(a), " genes)")
  centered <- similarity(ca, cb); scaled <- similarity(za, zb)
  half <- seq_len(nrow(a)) %% 2L == 0L
  split_a <- similarity(za[half, , drop = FALSE], zb[half, , drop = FALSE])
  split_b <- similarity(za[!half, , drop = FALSE], zb[!half, , drop = FALSE])
  rank_a <- apply(a, 2L, rank, ties.method = "average")
  rank_b <- apply(b, 2L, rank, ties.method = "average")
  spearman <- similarity(rank_a, rank_b)
  score <- scaled
  if (i == j) diag(score) <- -Inf
  nearest_a <- max.col(score, ties.method = "first")
  nearest_b <- max.col(t(score), ties.method = "first")
  idx <- if (i == j) which(upper.tri(score), arr.ind = TRUE) else
    as.matrix(expand.grid(seq_len(ncol(a)), seq_len(ncol(b))))
  aa <- idx[, 1]; bb <- idx[, 2]
  pairs <- data.table(study_a = sa, gsm_a = colnames(a)[aa],
    study_b = sb, gsm_b = colnames(b)[bb], genes = nrow(a),
    centered_r = centered[idx], standardized_r = scaled[idx],
    split_half_1_r = split_a[idx], split_half_2_r = split_b[idx],
    raw_spearman = spearman[idx],
    reciprocal_nearest = nearest_a[aa] == bb & nearest_b[bb] == aa)
  pairs[, same_gsm := gsm_a == gsm_b]
  pairs[, documented_relation := vapply(seq_len(.N), function(q)
    paste(sort(c(gsm_a[q], gsm_b[q])), collapse = "|") %in% relation_keys, TRUE)]
  # Broad retrieval threshold, not a validated duplicate decision threshold.
  pairs[, screen_flag := standardized_r >= cfg$screen_r &
    split_half_1_r >= cfg$screen_r & split_half_2_r >= cfg$screen_r]
  # Also retain each sample's nearest candidate to support manual review.
  pairs[, nearest_either := nearest_a[aa] == bb | nearest_b[bb] == aa]
  native_same <- registry$namespace[i] == registry$namespace[j]
  pairs[, `:=`(native_features = NA_integer_, native_equal_fraction = NA_real_,
                native_rmse = NA_real_, native_max_abs_diff = NA_real_)]
  if (native_same) {
    features <- intersect(rownames(da$native), rownames(db$native))
    selected <- which(pairs$screen_flag | pairs$same_gsm | pairs$documented_relation |
                        pairs$nearest_either)
    for (q in selected) {
      x <- da$native[features, pairs$gsm_a[q]]
      y <- db$native[features, pairs$gsm_b[q]]
      ok <- is.finite(x) & is.finite(y); diff <- x[ok] - y[ok]
      if (length(diff)) {
        set(pairs, q, "native_features", length(diff))
        set(pairs, q, "native_equal_fraction", mean(diff == 0))
        set(pairs, q, "native_rmse", sqrt(mean(diff^2)))
        set(pairs, q, "native_max_abs_diff", max(abs(diff)))
      }
    }
  }
  k <- k + 1L
  summaries[[k]] <- pairs[, .(study_a = sa, study_b = sb, genes = nrow(a),
    pairs_compared = .N, median_r = median(standardized_r),
    q99_r = as.numeric(quantile(standardized_r, .99)), max_r = max(standardized_r),
    screen_pairs = sum(screen_flag), known_same_gsm = sum(same_gsm),
    documented_pairs = sum(documented_relation))]
  correlation_rows[[k]] <- pairs[screen_flag | same_gsm | documented_relation | nearest_either]
}
hits <- rbindlist(correlation_rows)
hits <- merge(hits, samples[, .(study_a = study, gsm_a = gsm, title_a = title,
                               source_a = source)], by = c("study_a", "gsm_a"))
hits <- merge(hits, samples[, .(study_b = study, gsm_b = gsm, title_b = title,
                               source_b = source)], by = c("study_b", "gsm_b"))
setorder(hits, -standardized_r)
known <- hits[same_gsm | documented_relation]
new <- hits[!same_gsm & !documented_relation & screen_flag]
write_xlsx(list(inventory = rbindlist(inventory), pair_summary = rbindlist(summaries),
                new_candidates = new, known_overlaps = known,
                nearest_and_screened = hits),
           file.path(cfg$out_dir, "expression_overlap_audit.xlsx"))
capture.output(sessionInfo(), file = file.path(cfg$out_dir, "overlap_sessionInfo.txt"))
cat("\nNew high-similarity candidate pairs (not confirmed duplicates):", nrow(new), "\n")
print(head(new[, .(study_a, gsm_a, study_b, gsm_b, standardized_r, centered_r,
                   native_equal_fraction, native_rmse)], 40L))
cat("\nKnown overlaps recovered / above screen:\n")
print(known[, .(pairs = .N, screened = sum(screen_flag),
               min_r = min(standardized_r), max_r = max(standardized_r)),
            by = .(same_gsm, documented_relation)])
