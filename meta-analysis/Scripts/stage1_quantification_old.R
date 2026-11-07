#!/usr/bin/env Rscript
# =============================================================================
# MetaAnalysis — Stage 1: Raw data processing & uniform quantification
#
#   Affymetrix : CEL files -> RMA (RMA background correction, quantile
#                normalisation, median-polish summarisation, log2) with `oligo`,
#                one run per chip type, probe sets annotated to Entrez/Symbol.
#   RNA-Seq    : SRA FASTQ -> GRCh38 (GENCODE primary assembly) -> quantification
#                engine = "salmon"   : decoy-aware selective alignment to the
#                                      GRCh38 transcriptome, transcript-level
#                                      quantification, tximport to genes
#                engine = "rsubread" : splice-aware genome alignment to GRCh38
#                                      (Rsubread::align) + featureCounts
#                engine = "ncbi"     : NCBI GEO's uniformly re-processed counts
#                                      (HISAT2 -> GRCh38.p13 -> featureCounts)
#
# Why the series matrix is only the *manifest*:
#   A GEO series matrix holds values the submitters already processed (RMA with
#   their own CDF, quantile-normalised or voom counts, hg19, ...). Stage 1 must
#   start from raw data, so the series matrix is used for GSM IDs, CEL links,
#   SRX links and phenotype only. Its expression table is ignored.
#
# Layout expected (one folder per GSE):
#   <project>/Data/Raw/<GSE>/<GSE>[-GPLxxxx]_series_matrix.txt[.gz]
# Outputs:
#   <project>/Data/Stage1/<cohort>/...      final matrices (.rds + .tsv.gz), QC
#   <project>/Data/Stage1_work/<cohort>/... CEL / FASTQ / salmon / BAM files
#   <project>/Reference/GRCh38_gencode<rel>/ reference + indices (built once)
#
# Run on the Mac (needs access to ftp.ncbi.nlm.nih.gov and ftp.sra.ebi.ac.uk):
#   Rscript stage1_quantification.R                 # all cohorts
#   Rscript stage1_quantification.R GSE30122 GSE142025
# =============================================================================

# ---------------------------------------------------------------- 0. Config --
cfg <- list(
  project_dir   = path.expand("~/Striki/Quant/quant.2/MetaAnalysis"),
  threads       = max(1L, parallel::detectCores() - 2L),
  gencode_rel   = "46",        # GENCODE release to cite in Methods
  rnaseq_engine = "salmon",    # "salmon" | "rsubread" | "ncbi"
  keep_fastq    = FALSE,       # delete FASTQ after a sample is quantified
  rsubread_mem  = 8000,        # MB for Rsubread::buildindex
  qc            = TRUE         # raw/normalised boxplots + RLE/NUSE (arrays)
)
cfg$raw_dir  <- file.path(cfg$project_dir, "Data", "Raw")
cfg$out_dir  <- file.path(cfg$project_dir, "Data", "Stage1")
cfg$work_dir <- file.path(cfg$project_dir, "Data", "Stage1_work")
cfg$ref_dir  <- file.path(cfg$project_dir, "Reference",
                          paste0("GRCh38_gencode", cfg$gencode_rel))

args <- commandArgs(trailingOnly = TRUE)
cfg$cohorts <- if (length(args)) args else NULL   # NULL = everything in Data/Raw

options(timeout = 6 * 3600, stringsAsFactors = FALSE)
set.seed(1)

# ------------------------------------------------------------- 1. Packages --
need_cran <- c("BiocManager", "data.table")
need_bioc <- c("GEOquery", "Biobase", "oligo", "affyio", "AnnotationDbi",
               "tximport", "rtracklayer", "Rsubread")
# Platform design (pd.*) + annotation (*.db) packages for the chips in this project
need_affy <- c("pd.hg.u133a", "pd.hg.u133a.2", "pd.hg.u133.plus.2",
               "pd.hugene.2.1.st", "pd.hta.2.0",
               "hgu133a.db", "hgu133a2.db", "hgu133plus2.db",
               "hugene21sttranscriptcluster.db", "hta20transcriptcluster.db")

install_missing <- function() {
  miss <- setdiff(need_cran, rownames(installed.packages()))
  if (length(miss)) install.packages(miss)
  miss <- setdiff(c(need_bioc, need_affy), rownames(installed.packages()))
  if (length(miss)) BiocManager::install(miss, update = FALSE, ask = FALSE)
}
install_missing()
suppressPackageStartupMessages({
  library(data.table); library(Biobase); library(GEOquery); library(oligo)
})

# pd package (what oligo loads)  ->  annotation package (probe -> gene)
PD2DB <- c(
  "pd.hg.u133a"         = "hgu133a.db",
  "pd.hg.u133a.2"       = "hgu133a2.db",
  "pd.hg.u133.plus.2"   = "hgu133plus2.db",
  "pd.hugene.1.0.st.v1" = "hugene10sttranscriptcluster.db",
  "pd.hugene.2.0.st"    = "hugene20sttranscriptcluster.db",
  "pd.hugene.2.1.st"    = "hugene21sttranscriptcluster.db",
  "pd.hta.2.0"          = "hta20transcriptcluster.db"
)

msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)
dir_ok <- function(...) { d <- file.path(...); dir.create(d, TRUE, FALSE); d }
https <- function(u) sub("^ftp://", "https://", u)

run_cmd <- function(cmd, args) {
  msg(cmd, " ", paste(args, collapse = " "))
  st <- system2(cmd, args)
  if (st != 0) stop(cmd, " exited with status ", st)
}

fetch <- function(url, dest, md5 = NULL) {
  if (file.exists(dest) && (is.null(md5) || unname(tools::md5sum(dest)) == md5))
    return(invisible(dest))
  msg("download ", basename(dest))
  tmp <- paste0(dest, ".part")
  download.file(url, tmp, mode = "wb", quiet = TRUE)
  if (!is.null(md5) && unname(tools::md5sum(tmp)) != md5)
    stop("MD5 mismatch: ", basename(dest))
  file.rename(tmp, dest)
  invisible(dest)
}

# ------------------------------------------- 2. Series matrix -> manifest --
read_manifest <- function(sm_file) {
  es <- suppressMessages(GEOquery::getGEO(filename = sm_file, getGPL = FALSE))
  pd <- Biobase::pData(es)
  pd[] <- lapply(pd, as.character)
  pd
}

cols_like <- function(pd, pattern)
  unlist(pd[, grep(pattern, names(pd)), drop = FALSE], use.names = FALSE)

classify_cohort <- function(pd) {
  if (any(grepl("RNA-Seq", cols_like(pd, "^library_strategy"), ignore.case = TRUE)))
    return("rnaseq")
  if (any(grepl("\\.CEL(\\.gz)?$", cols_like(pd, "^supplementary_file"), ignore.case = TRUE)))
    return("affymetrix")
  "other"   # Illumina BeadChip, Agilent, ... -> not handled in this stage
}

# =============================================================================
# 3. AFFYMETRIX: CEL -> RMA
# =============================================================================
download_cels <- function(pd, dest) {
  supp <- pd[, grep("^supplementary_file", names(pd)), drop = FALSE]
  url  <- apply(supp, 1, function(x)
    x[grepl("\\.CEL(\\.gz)?$", x, ignore.case = TRUE)][1])
  if (anyNA(url)) stop("No CEL link for: ", paste(rownames(pd)[is.na(url)], collapse = ", "))
  files <- file.path(dest, basename(url))
  # Per-GSM download = only the samples in *this* series matrix. This also
  # handles multi-platform series (GSE99339, GSE104948, GSE104954) correctly.
  for (i in seq_along(url)) fetch(https(url[i]), files[i])
  setNames(files, rownames(pd))
}

annotate_probes <- function(ids, pd_pkg) {
  ann <- data.frame(PROBEID = ids, ENTREZID = NA_character_, SYMBOL = NA_character_)
  db  <- unname(PD2DB[pd_pkg])
  if (is.na(db)) { warning("No annotation package mapped for ", pd_pkg); return(ann) }
  if (!requireNamespace(db, quietly = TRUE)) BiocManager::install(db, update = FALSE, ask = FALSE)
  dbo <- getExportedValue(db, db)
  map <- function(col) tryCatch(
    unname(AnnotationDbi::mapIds(dbo, ids, col, "PROBEID", multiVals = "first")),
    error = function(e) rep(NA_character_, length(ids)))
  ann$ENTREZID <- map("ENTREZID")
  ann$SYMBOL   <- map("SYMBOL")
  ann
}

array_qc <- function(raw, eset, plm_args, pdf_file) {
  pdf(pdf_file, width = 12, height = 6)
  on.exit(dev.off())
  op <- par(mar = c(9, 4, 3, 1))
  oligo::boxplot(raw, which = "pm", las = 2, main = "Raw PM intensities (log2)")
  boxplot(exprs(eset), las = 2, outline = FALSE, main = "After RMA (log2)")
  plm <- do.call(oligo::fitProbeLevelModel, c(list(object = raw), plm_args))
  oligo::RLE(plm,  las = 2, main = "RLE")
  oligo::NUSE(plm, las = 2, main = "NUSE")
  par(op)
}

process_affymetrix <- function(cohort, pd, out, work) {
  cels <- download_cels(pd, dir_ok(work, "CEL"))

  # Group by chip type read from the CEL header: oligo cannot mix chips.
  chip <- vapply(cels, function(f) affyio::read.celfile.header(f)[["cdfName"]], "")
  res  <- list()

  for (ct in unique(chip)) {
    files <- cels[chip == ct]
    msg(cohort, ": RMA on ", length(files), " x ", ct)
    raw <- oligo::read.celfiles(filenames = files, sampleNames = names(files),
                                verbose = FALSE)
    is_st <- is(raw, "GeneFeatureSet") || is(raw, "ExonFeatureSet")

    # RMA = RMA convolution background correction + quantile normalisation
    #       + median-polish summarisation, output log2.
    # ST arrays (HuGene, HTA): summarise to transcript clusters ("core").
    rma_args <- list(object = raw, background = TRUE, normalize = TRUE)
    if (is_st) rma_args$target <- "core"
    eset <- do.call(oligo::rma, rma_args)

    # Drop Affymetrix control probe sets
    eset <- eset[!grepl("^AFFX", featureNames(eset), ignore.case = TRUE), ]

    pd_pkg <- annotation(raw)
    fData(eset) <- annotate_probes(featureNames(eset), pd_pkg)
    rownames(fData(eset)) <- featureNames(eset)
    pData(eset) <- pd[sampleNames(eset), , drop = FALSE]
    annotation(eset) <- pd_pkg
    experimentData(eset)@other$stage1 <- list(
      method = "oligo::rma", background = TRUE, normalize = "quantile",
      summary = "median polish", target = if (is_st) "core" else "probeset",
      chip = ct, pd = pd_pkg,
      pd_version = as.character(packageVersion(pd_pkg)))

    tag <- if (length(unique(chip)) > 1) paste0(cohort, "_", ct) else cohort
    saveRDS(eset, file.path(out, paste0(tag, "_RMA_eset.rds")))
    fwrite(data.table(fData(eset), exprs(eset)),
           file.path(out, paste0(tag, "_RMA_log2.tsv.gz")), sep = "\t")

    if (cfg$qc) tryCatch(
      array_qc(raw, eset, if (is_st) list(target = "core") else list(),
               file.path(out, paste0(tag, "_RMA_QC.pdf"))),
      error = function(e) warning(cohort, " QC failed: ", conditionMessage(e)))

    res[[ct]] <- data.frame(features = nrow(eset), samples = ncol(eset),
                            platform = pd_pkg)
    rm(raw); gc()
  }
  do.call(rbind, res)
}

# =============================================================================
# 4. RNA-SEQ: FASTQ -> GRCh38 -> quantification
# =============================================================================
srx_of <- function(pd) {
  rel <- pd[, grep("^relation", names(pd)), drop = FALSE]
  apply(rel, 1, function(x) {
    m <- regmatches(x, regexpr("SRX[0-9]+", x))
    if (length(m)) m[1] else NA_character_
  })
}

# ENA mirrors SRA and returns direct FASTQ links + MD5
ena_runs <- function(srx) {
  rbindlist(lapply(srx, function(x) {
    u <- paste0("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=", x,
                "&result=read_run&format=tsv&fields=",
                "run_accession,experiment_accession,library_layout,fastq_ftp,fastq_md5")
    fread(u, sep = "\t", colClasses = "character")
  }))
}

download_fastq <- function(run, dest) {
  url <- strsplit(run$fastq_ftp, ";")[[1]]
  md5 <- strsplit(run$fastq_md5, ";")[[1]]
  keep <- if (run$library_layout == "PAIRED") grepl("_[12]\\.fastq\\.gz$", url)
          else rep(TRUE, length(url))
  files <- file.path(dest, basename(url[keep]))
  for (i in seq_along(files)) fetch(paste0("https://", url[keep][i]), files[i], md5[keep][i])
  files
}

# ---- reference (built once, shared by all cohorts) ----
prepare_reference <- function(engine) {
  rd   <- dir_ok(cfg$ref_dir)
  rel  <- cfg$gencode_rel
  base <- sprintf("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_%s/", rel)
  f <- c(tx     = sprintf("gencode.v%s.transcripts.fa.gz", rel),
         genome = "GRCh38.primary_assembly.genome.fa.gz",
         gtf    = sprintf("gencode.v%s.primary_assembly.annotation.gtf.gz", rel))
  p <- setNames(file.path(rd, f), names(f))
  for (k in names(f)) fetch(paste0(base, f[[k]]), p[[k]])

  # transcript -> gene table from the same GTF
  t2g_file <- file.path(rd, "tx2gene.tsv")
  if (!file.exists(t2g_file)) {
    g <- rtracklayer::import(p[["gtf"]], feature.type = "transcript")
    fwrite(unique(data.table(TXNAME = g$transcript_id, GENEID = g$gene_id,
                             SYMBOL = g$gene_name, GENETYPE = g$gene_type)),
           t2g_file, sep = "\t")
  }
  ref <- c(as.list(p), tx2gene = t2g_file)

  if (engine == "salmon") {
    idx <- file.path(rd, "salmon_index")
    if (!file.exists(file.path(idx, "versionInfo.json"))) {
      # Decoy-aware index: whole GRCh38 genome appended as decoys, so reads from
      # unannotated loci are not forced onto transcripts.
      decoys <- file.path(rd, "decoys.txt")
      system(sprintf("gzip -dc '%s' | grep '^>' | cut -d ' ' -f1 | sed 's/^>//' > '%s'",
                     p[["genome"]], decoys))
      gentrome <- file.path(rd, "gentrome.fa.gz")
      system(sprintf("cat '%s' '%s' > '%s'", p[["tx"]], p[["genome"]], gentrome))
      run_cmd("salmon", c("index", "-t", gentrome, "-d", decoys, "-i", idx,
                          "--gencode", "-k", "31", "-p", cfg$threads))
    }
    ref$index <- idx
  }

  if (engine == "rsubread") {
    idx <- file.path(dir_ok(rd, "rsubread"), "GRCh38")
    if (!length(list.files(dirname(idx), "\\.reads$")))
      Rsubread::buildindex(basename = idx, reference = p[["genome"]],
                           gappedIndex = TRUE, indexSplit = TRUE,
                           memory = cfg$rsubread_mem)
    ref$index <- idx
  }
  ref
}

# ---- engine: salmon ----
quant_salmon <- function(fq, layout, out_dir, ref) {
  if (file.exists(file.path(out_dir, "quant.sf"))) return(file.path(out_dir, "quant.sf"))
  a <- c("quant", "-i", ref$index, "-l", "A", "-p", cfg$threads,
         "--seqBias", "-o", out_dir)
  if (layout == "PAIRED") {
    a <- c(a, "--gcBias",
           "-1", sort(fq[grepl("_1\\.fastq", fq)]),
           "-2", sort(fq[grepl("_2\\.fastq", fq)]))
  } else a <- c(a, "-r", fq)
  run_cmd("salmon", a)
  file.path(out_dir, "quant.sf")
}

summarise_salmon <- function(quant_files, ref, pd, cohort, out) {
  t2g <- fread(ref$tx2gene)
  # Gene-level: raw estimated counts + average tx length (for DESeq2/edgeR offsets)
  txi <- tximport::tximport(quant_files, type = "salmon",
                            tx2gene = t2g[, .(TXNAME, GENEID)],
                            countsFromAbundance = "no")
  # Gene-level: length-scaled TPM counts (for limma-voom without offsets)
  txi_ls <- tximport::tximport(quant_files, type = "salmon",
                               tx2gene = t2g[, .(TXNAME, GENEID)],
                               countsFromAbundance = "lengthScaledTPM")
  # Transcript-level
  txi_tx <- tximport::tximport(quant_files, type = "salmon", txOut = TRUE)

  gene_ann <- unique(t2g[, .(GENEID, SYMBOL, GENETYPE)])
  saveRDS(list(gene = txi, gene_lengthScaledTPM = txi_ls, transcript = txi_tx,
               gene_annotation = gene_ann, pheno = pd[colnames(txi$counts), , drop = FALSE]),
          file.path(out, paste0(cohort, "_salmon_tximport.rds")))
  ann <- gene_ann[match(rownames(txi$counts), GENEID)]
  fwrite(data.table(ann, round(txi$counts)),
         file.path(out, paste0(cohort, "_gene_counts.tsv.gz")), sep = "\t")
  fwrite(data.table(ann, txi$abundance),
         file.path(out, paste0(cohort, "_gene_TPM.tsv.gz")), sep = "\t")
  fwrite(data.table(TXNAME = rownames(txi_tx$abundance), txi_tx$abundance),
         file.path(out, paste0(cohort, "_transcript_TPM.tsv.gz")), sep = "\t")
  data.frame(features = nrow(txi$counts), samples = ncol(txi$counts),
             platform = paste0("salmon/GENCODE v", cfg$gencode_rel))
}

# ---- engine: rsubread (genome alignment + featureCounts) ----
align_rsubread <- function(fq, layout, bam, ref) {
  if (file.exists(bam)) return(bam)
  if (layout == "PAIRED") {
    r1 <- sort(fq[grepl("_1\\.fastq", fq)]); r2 <- sort(fq[grepl("_2\\.fastq", fq)])
  } else { r1 <- fq; r2 <- NULL }
  # Several runs per sample: align each, then merge by summing counts later.
  bams <- character()
  for (i in seq_along(r1)) {
    b <- if (length(r1) == 1) bam else sub("\\.bam$", sprintf("_run%d.bam", i), bam)
    Rsubread::align(index = ref$index, readfile1 = r1[i],
                    readfile2 = if (length(r2)) r2[i] else NULL,
                    type = "rna", output_file = b, nthreads = cfg$threads,
                    sortReadsByCoordinates = TRUE)
    bams <- c(bams, b)
  }
  bams
}

summarise_rsubread <- function(bams_by_gsm, paired, ref, pd, cohort, out) {
  all_bams <- unlist(bams_by_gsm)
  fc <- Rsubread::featureCounts(all_bams, annot.ext = ref$gtf,
                                isGTFAnnotationFile = TRUE,
                                GTF.featureType = "exon", GTF.attrType = "gene_id",
                                isPairedEnd = paired, countReadPairs = paired,
                                strandSpecific = 0,   # set 1/2 if the library is stranded
                                nthreads = cfg$threads)
  owner  <- rep(names(bams_by_gsm), lengths(bams_by_gsm))
  counts <- sapply(split(seq_along(owner), owner),
                   function(j) rowSums(fc$counts[, j, drop = FALSE]))
  t2g <- unique(fread(ref$tx2gene)[, .(GENEID, SYMBOL, GENETYPE)])
  saveRDS(list(counts = counts, annotation = fc$annotation, stat = fc$stat,
               pheno = pd[colnames(counts), , drop = FALSE]),
          file.path(out, paste0(cohort, "_rsubread_featureCounts.rds")))
  fwrite(data.table(t2g[match(rownames(counts), GENEID)], counts),
         file.path(out, paste0(cohort, "_gene_counts.tsv.gz")), sep = "\t")
  data.frame(features = nrow(counts), samples = ncol(counts),
             platform = paste0("Rsubread/GENCODE v", cfg$gencode_rel))
}

# ---- engine: NCBI GEO uniformly processed counts ----
process_ncbi_counts <- function(gse, pd, cohort, out) {
  base <- "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&format=file"
  cf <- sprintf("%s_raw_counts_GRCh38.p13_NCBI.tsv.gz", gse)
  af <- "Human.GRCh38.p13.annot.tsv.gz"
  ok <- tryCatch({
    fetch(sprintf("%s&acc=%s&file=%s", base, gse, cf), file.path(out, cf))
    fetch(sprintf("%s&file=%s", base, af), file.path(out, af)); TRUE
  }, error = function(e) FALSE)
  if (!ok) stop("NCBI has no GRCh38 counts for ", gse,
                " (no public SRA reads?) - use engine 'salmon' or 'rsubread'")
  cnt  <- fread(file.path(out, cf))
  keep <- intersect(rownames(pd), names(cnt))      # this platform's GSMs only
  cnt  <- cnt[, c("GeneID", keep), with = FALSE]
  ann  <- fread(file.path(out, af))
  saveRDS(list(counts = as.matrix(cnt[, -1], rownames = cnt$GeneID),
               annotation = ann[match(cnt$GeneID, GeneID)],
               pheno = pd[keep, , drop = FALSE]),
          file.path(out, paste0(cohort, "_NCBI_GRCh38_counts.rds")))
  data.frame(features = nrow(cnt), samples = length(keep),
             platform = "NCBI HISAT2/featureCounts GRCh38.p13")
}

process_rnaseq <- function(cohort, gse, pd, out, work) {
  engine <- cfg$rnaseq_engine
  if (engine == "ncbi") return(process_ncbi_counts(gse, pd, cohort, out))

  srx <- srx_of(pd)
  if (all(is.na(srx)))
    stop("No public SRA reads linked (e.g. controlled access) - raw-read ",
         "re-quantification impossible; keep the deposited matrix and flag it.")
  if (anyNA(srx)) warning(cohort, ": no SRX for ", paste(names(srx)[is.na(srx)], collapse = ", "))
  srx  <- srx[!is.na(srx)]
  runs <- ena_runs(unique(srx))
  fwrite(runs, file.path(out, paste0(cohort, "_SRA_runs.tsv")), sep = "\t")

  ref    <- prepare_reference(engine)
  fq_dir <- dir_ok(work, "FASTQ")
  res    <- list()

  for (gsm in names(srx)) {
    r      <- runs[experiment_accession == srx[[gsm]]]
    layout <- unique(r$library_layout)
    if (length(layout) != 1) stop(gsm, ": mixed library layouts")
    target <- if (engine == "salmon") file.path(work, "salmon", gsm, "quant.sf")
              else file.path(work, "bam", paste0(gsm, ".bam"))
    msg(cohort, " ", gsm, " (", nrow(r), " run(s), ", layout, ")")
    if (!file.exists(target)) {
      fq <- unlist(lapply(split(r, r$run_accession), download_fastq, dest = fq_dir))
      res[[gsm]] <- if (engine == "salmon")
        quant_salmon(fq, layout, file.path(dir_ok(work, "salmon"), gsm), ref)
      else align_rsubread(fq, layout, file.path(dir_ok(work, "bam"), paste0(gsm, ".bam")), ref)
      if (!cfg$keep_fastq) unlink(fq)
    } else {
      res[[gsm]] <- if (engine == "salmon") target
                    else list.files(dirname(target), paste0("^", gsm, "(_run[0-9]+)?\\.bam$"),
                                    full.names = TRUE)
    }
  }

  if (engine == "salmon")
    summarise_salmon(unlist(res), ref, pd, cohort, out)
  else
    summarise_rsubread(res, all(runs$library_layout == "PAIRED"), ref, pd, cohort, out)
}

# =============================================================================
# 5. Driver
# =============================================================================
sm_files <- list.files(cfg$raw_dir, "_series_matrix\\.txt(\\.gz)?$",
                       recursive = TRUE, full.names = TRUE)
cohorts  <- sub("_series_matrix\\.txt(\\.gz)?$", "", basename(sm_files))
sm_files <- setNames(sm_files, cohorts)[!duplicated(cohorts)]   # .txt preferred over .gz
if (!is.null(cfg$cohorts))
  sm_files <- sm_files[sub("-GPL.*$", "", names(sm_files)) %in% cfg$cohorts |
                       names(sm_files) %in% cfg$cohorts]

dir_ok(cfg$out_dir)
log <- list()

for (cohort in names(sm_files)) {
  gse  <- sub("-GPL.*$", "", cohort)
  pd   <- read_manifest(sm_files[[cohort]])
  type <- classify_cohort(pd)
  out  <- dir_ok(cfg$out_dir, cohort)
  work <- dir_ok(cfg$work_dir, cohort)
  msg("==== ", cohort, " [", type, ", ", nrow(pd), " samples]")
  fwrite(data.table(GSM = rownames(pd), pd),
         file.path(out, paste0(cohort, "_pheno.tsv")), sep = "\t")

  r <- tryCatch(
    switch(type,
           affymetrix = process_affymetrix(cohort, pd, out, work),
           rnaseq     = process_rnaseq(cohort, gse, pd, out, work),
           other      = stop("not Affymetrix or RNA-Seq - needs its own normalisation path")),
    error = function(e) { warning(cohort, ": ", conditionMessage(e)); conditionMessage(e) })

  log[[cohort]] <- if (is.data.frame(r))
    data.frame(cohort, type, status = "done", r)
  else
    data.frame(cohort, type, status = "skipped/failed", features = NA,
               samples = nrow(pd), platform = r)
}

summary_tab <- rbindlist(log, fill = TRUE)
fwrite(summary_tab, file.path(cfg$out_dir, "stage1_summary.tsv"), sep = "\t")
print(summary_tab)
writeLines(capture.output(sessionInfo()), file.path(cfg$out_dir, "stage1_sessionInfo.txt"))
if (cfg$rnaseq_engine == "salmon" && nzchar(Sys.which("salmon")))
  writeLines(system2("salmon", "--version", stdout = TRUE),
             file.path(cfg$out_dir, "stage1_salmon_version.txt"))
