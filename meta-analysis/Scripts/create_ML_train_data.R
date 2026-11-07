#!/usr/bin/env Rscript
# =============================================================================
# RenoRisk: sample-level ML input, before cross-study batch correction.
# Statistical preprocessing must be fitted inside training/validation folds.
# Microarray + RNA-Seq in one table (by instruction, 2026-11-03), `Platform`
# column marks each row. GSE142025 comes through Stage23 like every cohort;
# GSE137570 and GSE219185 (no controls) come from Scripts/by-claude/rnaseq_panel.R and
# enter the ML inputs only. Run after Stage23, Stage4 and rnaseq_panel.R.
# =============================================================================
set.seed(1)
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(writexl)
})

# 1. Configuration and source files.
cfg <- list(
  project_dir = path.expand("~/Striki/Quant/quant.2/MetaAnalysis"),
  genes = c("C4BPA", "COL13A1", "HABP2", "CFH", "CFI"),
  excluded_studies = c("GSE30528", "GSE30529"),
  rnaseq = c("GSE142025", "GSE141295", "GSE175759", "GSE137570-GPL23227", "GSE219185")   # Platform = "RNA-Seq"
)
cfg$rnaseq_dir <- file.path(cfg$project_dir, "Results", "RNASeq")
cfg$source_dir <- file.path(cfg$project_dir, "Results", "Stage23")
cfg$sex_file <- file.path(cfg$project_dir, "Results", "Stage4", "sex_by_sample.xlsx")
cfg$out_dir <- file.path(cfg$project_dir, "Results", "ML_training")
cfg$ml_input_dir <- file.path(cfg$project_dir, "Data", "ML_input")   # all.xlsx: ML input
for (d in c(cfg$out_dir, cfg$ml_input_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

expression <- as.data.table(read_excel(
  file.path(cfg$source_dir, "panel_uncorrected_long.xlsx"),
  col_types = c(rep("text", 6L), "numeric")))
indications <- as.data.table(read_excel(
  file.path(cfg$source_dir, "sample_indication.xlsx"),
  col_types = rep("text", 6L)))
sex_data <- as.data.table(read_excel(cfg$sex_file, col_types = rep("text", 6L)))

# RNA-Seq series without controls (rnaseq_panel.R): ML inputs only.
rna_long <- as.data.table(read_excel(file.path(cfg$rnaseq_dir, "panel_rnaseq_long.xlsx"),
                                     col_types = c(rep("text", 6L), "numeric")))
rna_smp  <- as.data.table(read_excel(file.path(cfg$rnaseq_dir, "rnaseq_samples.xlsx"),
                                     col_types = rep("text", 7L)))
expression  <- rbind(expression, rna_long)
indications <- rbind(indications, rna_smp[, .(StudyID, GSM, Compartment, Condition,
                                              Indication, Label)])
sex_data <- rbind(sex_data, rna_smp[, .(StudyID, GSM, Sex,
  SexSource = fifelse(is.na(Sex), "unknown", "rnaseq_panel"))], fill = TRUE)

# 2. Exclude the two subset series, retaining their samples through GSE30122.
#    Keep other suspected overlaps for review; GSM uniqueness is not proof of
#    independent patients. Do not alter existing accepted cohort inclusion rules.
expression <- expression[!StudyID %in% cfg$excluded_studies]
indications <- indications[!StudyID %in% cfg$excluded_studies]
sex_data <- sex_data[!StudyID %in% cfg$excluded_studies]
stopifnot(
  setequal(unique(expression$Gene), cfg$genes),
  all(is.finite(expression$Expression)),
  all(expression$Condition %in% c("Control", "Disease")),
  !anyDuplicated(expression[, .(StudyID, GSM, Gene)]),
  !anyDuplicated(indications[, .(StudyID, GSM)]),
  !anyDuplicated(sex_data[, .(StudyID, GSM)]),
  all(is.na(sex_data$Sex) | sex_data$Sex %in% c("M", "F")),
  all(sex_data$SexSource %in% c("reported", "inferred", "unknown", "rnaseq_panel"))
)

# 3. One row per GSM. Preserve numeric precision and genuine missing metadata.
training_data <- dcast(
  expression, GSM + StudyID + Batch + Compartment + Condition ~ Gene,
  value.var = "Expression")
training_data <- merge(
  training_data, indications[, .(StudyID, GSM, Indication)],
  by = c("StudyID", "GSM"), all.x = TRUE, sort = FALSE)
# Prefer the Stage4 reported value when present; otherwise retain its existing
# expression-based inference. Sex stays as Stage4 codes it: "M" / "F" / NA.
stopifnot(nrow(training_data[!sex_data, on = .(StudyID, GSM)]) == 0L)
training_data <- merge(
  training_data, sex_data[, .(StudyID, GSM, Sex)],
  by = c("StudyID", "GSM"), all.x = TRUE, sort = FALSE)
training_data[, Disease_Label := as.integer(Condition == "Disease")]
training_data[, Platform := fifelse(StudyID %in% cfg$rnaseq, "RNA-Seq", "Microarray")]
columns <- c("GSM", "StudyID", "Batch", "Compartment", "Condition",
             "Disease_Label", "Indication", "Sex", "Platform", cfg$genes)
setcolorder(training_data, columns)
setorder(training_data, StudyID, GSM)
# A gene is NA only where the study's deposit lacks it (C4BPA: GSE142025,
# GSE137570); every value present is finite.
present <- unique(expression[, .(StudyID, Gene)])
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

# 4. Keep the feature table first and explain label/provenance choices separately.
readme <- data.table(
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
    "GSE30528 and GSE30529 (all samples in GSE30122) and GSE35487 (same arrays as GSE35488) and GSE37460 (arrays already in GSE104948/GSE99339) are not used.",
    "Unique GSM accessions are not necessarily independent biological samples or patients. Documented reanalyses remain included.",
    "The expression audit flags GSE108109/GSE200828 as strongly overlapping and GSE121233/GSE129973 as additional candidates. See Results/SampleOverlap/README.md. GSE219185 (NEPTUNE tubulointerstitium) states it shares samples with NEPTUNE GEO sets including GSE108113, whose glomerular part is GSE108109: same patients, different compartment.",
    "Use study/patient-aware validation (leave-one-study-out). Keep linked studies together (e.g. GSE219185 with GSE108109); fit learned preprocessing only within training folds. Sex is coded M / F (categorical). No fixed train/test split is assigned here.",
    "Expression/indications: Results/Stage23/panel_uncorrected_long.xlsx and sample_indication.xlsx. Sex: Results/Stage4/sex_by_sample.xlsx. GSE137570/GSE219185: Results/RNASeq/ (rnaseq_panel.R). Submitter-processed values from Data/Raw_filtered only.",
    "Microarray (931 rows) or RNA-Seq (391 rows: GSE142025 36, GSE141295 24, GSE175759 87, GSE137570-GPL23227 24, GSE219185 220). Platform is confounded with study: use it to stratify, not as a predictor.",
    "C4BPA is NA for GSE142025 and GSE137570 (60 rows): not in those deposits. Impute inside training folds or drop C4BPA; most scikit-learn estimators reject NaN.",
    "GSE137570 (CKD) and GSE219185 (FSGS/MCD) are all Disease. Their study identity predicts the label, so only grouped-by-study CV is valid; random-split CV is optimistic."
  ))
output <- file.path(cfg$ml_input_dir, "all.xlsx")
write_xlsx(list(training_data = training_data, README = readme), output)

# 5. Verify that the written workbook preserves feature values and sample labels.
reloaded <- as.data.table(read_excel(output, sheet = "training_data",
  col_types = c(rep("text", 5L), "numeric", "text", "text", "text", rep("numeric", 5L))))
stopifnot(
  identical(reloaded$GSM, training_data$GSM),
  identical(as.matrix(reloaded[, cfg$genes, with = FALSE]),
            as.matrix(training_data[, cfg$genes, with = FALSE])),
  identical(as.integer(reloaded$Disease_Label), training_data$Disease_Label),
  identical(reloaded$Sex, training_data$Sex)
)
cat("Created", output, "\n")

# 6. Same rows and metadata, batch-corrected. The Stage23 model
#    (limma::removeBatchEffect; plan "~ StudyID + Condition": batch = StudyID,
#    design = intercept + Disease), refitted here over ALL ML rows so the
#    RNA-Seq series without controls (not in Stage23) are corrected too. Gene by
#    gene on the rows that have the gene. A disease-only study's offset rests on
#    the assumption that its disease effect equals the pooled one.
#    The labels shaped these values: visualisation, not unbiased validation.
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

readme_bn <- copy(readme)
readme_bn[Item == "Expression scale", Description :=
  "limma::removeBatchEffect with the Stage23 model (plan: ~ StudyID + Condition; batch = StudyID, 19 studies; design = intercept + Disease), refitted over all 1322 rows gene by gene, so values differ from Results/Stage23/panel_harmonised_long.xlsx (fitted on 1078 rows). Offsets of the two disease-only RNA-Seq series (GSE137570, GSE219185) assume their disease effect equals the pooled one. Cross-compartment level differences are removed except within GSE30122 (one study, both compartments). Gene-centred cohorts make the common scale relative, not absolute."]
readme_bn[Item == "Batch correction", Description :=
  "Corrected for batch location only, not spread: residual per-gene SD still differs up to ~6.7x between studies. The correction design protected a Disease effect, so Disease_Label influenced these values."]
readme_bn[Item == "Validation", Description :=
  "LABEL LEAKAGE: the correction was fitted on all samples with their disease labels, including any later test fold. Cross-validated performance on these values is optimistic. Use Data/ML_input/all.xlsx with in-fold preprocessing for unbiased validation; use this workbook for visualisation or comparison. Sex is coded M / F (categorical)."]
readme_bn[Item == "Input provenance", Description :=
  "Expression: Data/ML_input/all.xlsx values, corrected in create_ML_train_data.R section 6. Submitter-processed values from Data/Raw_filtered only."]
output_bn <- file.path(cfg$out_dir, "all_bn.xlsx")
write_xlsx(list(training_data = batch_normalized, README = readme_bn), output_bn)

reloaded_bn <- as.data.table(read_excel(output_bn, sheet = "training_data",
  col_types = c(rep("text", 5L), "numeric", "text", "text", "text", rep("numeric", 5L))))
stopifnot(
  identical(reloaded_bn$GSM, batch_normalized$GSM),
  isTRUE(all.equal(as.matrix(reloaded_bn[, cfg$genes, with = FALSE]),
                   as.matrix(batch_normalized[, cfg$genes, with = FALSE]), tolerance = 1e-12)),
  identical(as.integer(reloaded_bn$Disease_Label), batch_normalized$Disease_Label)
)
capture.output(sessionInfo(),
               file = file.path(cfg$out_dir, "ML_training_sessionInfo.txt"))
cat("Created", output_bn, "\n")
print(training_data[, .N, by = .(Condition)])
print(training_data[, .N, by = .(Compartment)])
print(training_data[, .N, by = .(Sex)])
