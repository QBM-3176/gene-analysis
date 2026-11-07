#!/usr/bin/env Rscript
# =============================================================================
# RenoRisk: all.xlsx with ComBat (sva) cross-study batch correction -> all_corr.xlsx.
# Same rows, columns and metadata as all.xlsx; only the five gene columns change.
# Plan §3 Stage 3, "~ StudyID + Condition" (by instruction, 2026-11-04): batch =
# StudyID (19 studies), mod = intercept + Disease, as Stage23 does. Compartment is
# aliased with study (except GSE30122). Labels shape the values, so CV on this
# file is optimistic. Two RNA-Seq series have no controls (GSE137570, GSE219185):
# their offsets assume their disease effect equals the pooled one. C4BPA is absent from two RNA-Seq units, and ComBat cannot fit a
# batch with no values for a gene, so it runs in two passes: the complete genes
# over all rows; each incomplete gene with all genes over the rows that have it
# (only that gene is taken from the second pass; EB priors still use all genes).
# Usage : Rscript Scripts/create_ML_corr_data.R   (run after create_ML_train_data.R)
# =============================================================================
set.seed(1)
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(writexl)
})

cfg <- list(
  project_dir = path.expand("~/Striki/Quant/quant.2/MetaAnalysis"),
  genes = c("C4BPA", "COL13A1", "HABP2", "CFH", "CFI")
)
cfg$ml_input_dir <- file.path(cfg$project_dir, "Data", "ML_input")
cfg$out_dir <- file.path(cfg$project_dir, "Results", "ML_training")
col_types <- c(rep("text", 5L), "numeric", "text", "text", "text", rep("numeric", 5L))

input <- file.path(cfg$ml_input_dir, "all.xlsx")
data <- as.data.table(read_excel(input, sheet = "training_data", col_types = col_types))
readme <- as.data.table(read_excel(input, sheet = "README", col_types = c("text", "text")))
stopifnot(nrow(data) == 1322L, !anyDuplicated(data$GSM))

# ComBat: genes x samples, parametric empirical Bayes location + scale adjustment.
combat <- function(rows, genes) {
  d <- data[rows]
  x <- t(as.matrix(d[, genes, with = FALSE])); colnames(x) <- d$GSM
  sva::ComBat(x, batch = d$StudyID, mod = model.matrix(~ Disease_Label, data = d),   # plan §3: ~ StudyID + Condition
              par.prior = TRUE)
}
corr <- copy(data)
complete <- cfg$genes[vapply(cfg$genes, function(g) !anyNA(data[[g]]), TRUE)]
corr[, (complete) := as.data.table(t(combat(seq_len(nrow(data)), complete)))]
for (g in setdiff(cfg$genes, complete)) {
  rows <- which(!is.na(data[[g]]))
  stopifnot(!anyNA(data[rows, cfg$genes, with = FALSE]))
  set(corr, rows, g, as.numeric(combat(rows, cfg$genes)[g, ]))
}
stopifnot(identical(corr$GSM, data$GSM),
          all(vapply(cfg$genes, function(g)
            identical(is.na(corr[[g]]), is.na(data[[g]])) && all(is.finite(corr[[g]][!is.na(corr[[g]])])), TRUE)))

# Share of each gene's variance explained by the batch unit, before vs after.
batch_r2 <- function(d) sapply(cfg$genes, function(g)
  summary(lm(d[[g]] ~ d$Batch))$r.squared)   # lm drops the NA rows
r2 <- data.table(Gene = cfg$genes, R2_batch_before = round(batch_r2(data), 3),
                 R2_batch_after = round(batch_r2(corr), 3))
print(r2)

readme[Item == "Expression scale", Description :=
  "all.xlsx values after sva::ComBat (parametric empirical Bayes, location and scale) over all 1322 samples (931 microarray, 391 RNA-Seq), batch = StudyID (19 studies; plan: ~ StudyID + Condition); C4BPA in a second pass over the 1262 rows that have it. Cross-compartment level differences are removed except within GSE30122 (one study, both compartments). Gene-centred cohorts make the common scale relative, not absolute."]
readme[Item == "Batch correction", Description :=
  "ComBat adjusts batch mean and variance per gene. Its empirical Bayes priors are estimated from only 5 genes (4 in the first pass), so the shrinkage is weak. Offsets of the disease-only RNA-Seq series (GSE137570, GSE219185) assume their disease effect equals the pooled one. The design protected a Disease effect, so Disease_Label influenced these values."]
readme[Item == "Validation", Description :=
  "LABEL LEAKAGE: the correction was fitted on all samples with their disease labels, including any later test fold. Cross-validated performance on these values is optimistic. Use Data/ML_input/all.xlsx with in-fold preprocessing for unbiased validation. Sex is coded M / F (categorical)."]
readme[Item == "Input provenance", Description :=
  "Data/ML_input/all.xlsx (create_ML_train_data.R), corrected by Scripts/create_ML_corr_data.R. Submitter-processed values (microarray and RNA-Seq) from Data/Raw_filtered only."]

output <- file.path(cfg$ml_input_dir, "all_corr.xlsx")
write_xlsx(list(training_data = corr, README = readme, batch_r2 = r2), output)

reloaded <- as.data.table(read_excel(output, sheet = "training_data", col_types = col_types))
stopifnot(identical(reloaded$GSM, corr$GSM),
          isTRUE(all.equal(as.matrix(reloaded[, cfg$genes, with = FALSE]),
                           as.matrix(corr[, cfg$genes, with = FALSE]), tolerance = 1e-12)),
          identical(as.integer(reloaded$Disease_Label), as.integer(corr$Disease_Label)))
capture.output(sessionInfo(), file = file.path(cfg$out_dir, "all_corr_sessionInfo.txt"))
cat("Created", output, ":", nrow(corr), "samples\n")
