#!/usr/bin/env Rscript
# =============================================================================
# MetaAnalysis — Stage 4: multi-factor statistical modelling of the 5-gene panel
#
# SOURCE: CKD_Biomarker_Pooled_Analysis_Plan_Final.pdf — "PROJECT ANALYSIS PLAN &
#   TECHNICAL MASTER PLAN", RenoRisk Validation Across Multiple CKD Indications,
#   Version 2.0 (Final Analysis Plan), July 2026. Project Lead: Amela Jusic;
#   Execution Lead: Stojanko Tudjarski. Quotations below are verbatim.
#
# PLAN REQUIREMENTS IMPLEMENTED HERE
#   §3 Stage 4 — Multi-Factor Statistical Modeling:
#     "Model differential expression of the 5 RenoRisk genes controlling for
#     tissue compartment (glomerular vs. tubulointerstitial) and patient clinical
#     covariates."
#   §3 Stage 3 (model the plan names for the study/condition structure):
#     "... modeling ~ StudyID + Condition."
#   §1 Project Roles (Project Lead responsibilities):
#     "Translational alignment with clinical endpoints (eGFR, %TIF, CKD stages)"
#   §2 Target RenoRisk Biomarker Panel Specification:
#     "The RenoRisk panel consists of 5 specific target biomarkers categorized into
#     primary diagnostic candidates and additional complement system regulators"
#     Main RenoRisk: C4BPA (C4b-Binding Protein Alpha), COL13A1 (Collagen Type XIII
#     Alpha 1 Chain), HABP2 (Hyaluronan Binding Protein 2); Additional RenoRisk:
#     CFH (Complement Factor H), CFI (Complement Factor I).
#   Annex 2, column "Confirmed Quantitative Clinical Variables" (as stated there):
#     GSE66494  "eGFR and %TIF (tubulointerstitial fibrosis) - study built around
#               these two outcomes"
#     GSE104954 "eGFR (CKD-EPI) confirmed present in per-sample characteristics
#               field; diagnosis; technical outlier flag"
#     GSE104948 "eGFR, diagnosis - companion glomerular series to GSE104954, same
#               patients"
#     GSE129973 "Diagnosis-level metadata confirmed; per-sample proteinuria/eGFR
#               NOT confirmed in GEO record itself"
#     GSE121233 "Small n limits quantitative modeling"
#     GSE108109 "Full NEPTUNE clinical set (eGFR, UPCR, IFTA score) available, but
#               mainly via Nephroseq/NEPTUNE portal, not fully embedded in GEO"
#     GSE99339  "MEST-C histologic score correlatable per publications"
#     others    "Not confirmed in metadata - verify directly on GEO"
#
# HOW THIS SCRIPT ADDRESSES THEM (and where it deviates, with the reason)
#   * One weighted model per gene over the 967 aggregated samples (GSE142025 RNA-Seq included):
#       Expression ~ 0 + Batch + Disease:Compartment + Sex
#     StudyID enters as the batch unit (study x compartment, as in Stage 3);
#     Condition enters as a Disease effect per compartment.
#   * "controlling for tissue compartment": compartment is nested in batch, so its
#     main effect is aliased with the batch terms. The script therefore reports
#     the Disease effect within each compartment (Glomerular, Tubulointerstitial,
#     and the third stratum Whole biopsy, which the plan does not name) and the
#     between-compartment differences of that effect.
#   * "patient clinical covariates": every series matrix was checked. None of the
#     13 carries eGFR, %TIF, proteinuria, CKD stage or age — including ERCB and
#     GSE66494, for which Annex 2 lists eGFR/%TIF: those values live in
#     Nephroseq, not in GEO. Sex is the only available covariate: reported for
#     GSE116626, GSE129973, GSE96804; inferred elsewhere from XIST vs Y-chromosome
#     genes (validated against the reported cohorts, sex_by_cohort.xlsx).
#   * Input is the Stage 2/3 values BEFORE removeBatchEffect: limma advises
#     against modelling batch-corrected values; batch is modelled explicitly
#     instead. Weights = 1 / residual variance per batch unit, because spread
#     differs up to ~6.7x between studies.
#
# Plan: "Model differential expression of the 5 RenoRisk genes controlling for
# tissue compartment (glomerular vs. tubulointerstitial) and patient clinical
# covariates."
#
# Input : Results/Stage23/panel_uncorrected_long.xlsx  (written by Stage 2/3 —
#         same samples, labels, MaxMean collapse and log2 step, but BEFORE
#         removeBatchEffect) + the series matrices, for sex markers only (for the
#         RNA-Seq cohort GSE142025, its <id>_expression.txt).
# Output: Results/Stage4/ (all tables as .xlsx, the figure and sessionInfo)
#
# MODEL (one per gene, 967 samples)
#   Expression ~ 0 + Batch + Disease:Compartment + Sex
#   * Batch is fitted, not pre-removed: limma advises against running linear
#     models on removeBatchEffect output. Batch = study x compartment unit, as
#     in Stage 2/3 (GSE30122 split in two -> 12 units).
#   * Compartment main effects are NOT estimable — compartment is nested in the
#     batch units — so compartment is controlled by the batch terms, and the
#     Disease effect is estimated separately per compartment. Differences in
#     the disease effect BETWEEN compartments are estimable and tested.
#   * Weighted least squares: residual spread differs up to ~6.7x between
#     studies (Stage 2/3 residual_sd_ratio.xlsx), which breaks OLS's equal-
#     variance assumption. Two-step FGLS: OLS fit -> residual variance per
#     batch unit -> refit with weights 1/variance.
#
# CLINICAL COVARIATES — what the data allows
#   None of the 13 series matrices carries eGFR, %TIF, proteinuria or age.
#   Sex is the only patient covariate, and is REPORTED for 3 cohorts only
#   (GSE116626, GSE129973, GSE96804). For the rest it is INFERRED from
#   expression: XIST vs Y-linked genes (RPS4Y1, DDX3Y, KDM5D, UTY, EIF1AY).
#   The inference is validated against the reported cohorts (sex_by_cohort.xlsx).
#
# CAVEATS
#   * GSE30528/GSE30529 (cfg$subsets) are subsets of GSE30122. They are used in
#     the single-dataset calculations (sex per cohort, per-study effects) but
#     left out of the multi-factor model, which aggregates across datasets, so
#     no sample is counted twice.
#   * GSE30122 glomerular and tubular samples partly come from the same
#     individuals; that pairing is not modelled.
#   * Disease is the binary Control/Disease collapse of Stage 2/3.
#   * Values are submitter-processed (five different pipelines), and four
#     cohorts are per-gene centred — effects are relative within each study.
# =============================================================================

cfg <- list(
  project_dir = path.expand("~/Striki/Quant/quant.2/MetaAnalysis"),
  genes       = c("C4BPA", "COL13A1", "HABP2", "CFH", "CFI"),
# GSE30528/GSE30529 (subsets of GSE30122) were dropped and renamed -dup on 2026-11-03 (by instruction); the guard stays for any future contained dataset.
  subsets     = character(0),   # datasets fully contained in another: single-dataset use only
  y_genes     = c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY"),
  x_gene      = "XIST",
  min_sex_sep = 4,      # inferred sex used only if the two clusters are this many
                        # within-cluster SDs apart; otherwise left unknown
  width       = 9,
  height      = 5,
  dpi         = 300
)
cfg$results_dir <- file.path(cfg$project_dir, "Results")   # everything the pipeline generates
cfg$in_file <- file.path(cfg$results_dir, "Stage23", "panel_uncorrected_long.xlsx")
cfg$raw_dir <- file.path(cfg$project_dir, "Data", "Raw_filtered")
cfg$out_dir <- file.path(cfg$results_dir, "Stage4")
cfg$fig_dir <- cfg$out_dir

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(AnnotationDbi); library(WGCNA)
  library(readxl); library(writexl)
})
set.seed(1)
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

if (!file.exists(cfg$in_file))
  stop("missing ", cfg$in_file, " - run Scripts/stage23_panel_boxplots.R first")

# Probe-ID route per cohort (same as the Stage 2/3 registry), and Entrez IDs
# for the custom-CDF cohorts that address probes as "<entrez>_at".
IDTYPE <- c(`GSE104948-GPL22945` = "entrez", `GSE104954-GPL22945` = "entrez",
            GSE108109 = "entrez", GSE200828 = "entrez", `GSE99339-GPL19184` = "entrez",
            GSE121233 = "hta20transcriptcluster.db", GSE129973 = "hta20transcriptcluster.db",
            GSE96804  = "hta20transcriptcluster.db",
            GSE93798 = "entrez", GSE35488 = "entrez", GSE111154 = "hta20transcriptcluster.db",
            GSE30122 = "hgu133a2.db",
            GSE116626 = "illuminaHumanv4.db", GSE66494 = "hgug4112a.db",
            GSE142025 = "symbol", GSE141295 = "symbol", GSE175759 = "symbol")   # RNA-Seq: rows are gene symbols
# RNA-Seq: the series-matrix table is empty; values are in this file.
EXPR_FILE <- c(GSE142025 = "GSE142025_expression.txt",
               GSE141295 = "GSE141295_expression.txt",   # DESeq2 VST from counts
               GSE175759 = "GSE175759_expression.txt")   # (rnaseq_counts_vst.R)
LINEAR <- "GSE66494"   # Agilent FeatureExtraction output, needs log2
ENTREZ_SEX <- c(XIST = "7503", RPS4Y1 = "6192", DDX3Y = "8653", KDM5D = "8284",
                UTY = "7404", EIF1AY = "9086")

# ---------------------------------------------------------------------------
# 1. Panel values (uncorrected) from Stage 2/3
# ---------------------------------------------------------------------------
# Column types are fixed, not guessed: GSM, StudyID, Condition, Compartment,
# Batch, Gene (text), Expression (numeric).
long <- as.data.table(read_excel(cfg$in_file, col_types = c(rep("text", 6), "numeric")))
long[, Gene := factor(Gene, levels = intersect(cfg$genes, unique(Gene)))]
samples <- unique(long[, .(GSM, StudyID, Condition, Compartment, Batch)])
msg(nrow(samples), " observations (", uniqueN(samples$GSM), " distinct samples), ",
    uniqueN(samples$Batch), " batch units, ", nlevels(long$Gene), " genes")

# ---------------------------------------------------------------------------
# 2. Sex: reported where the cohort records it, else inferred from expression
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
  b <- which(trimws(lines) == "!series_matrix_table_begin") + 1
  e <- which(trimws(lines) == "!series_matrix_table_end")   - 1
  expr <- as.matrix(read.table(text = paste(lines[b:e], collapse = "\n"),
                               header = TRUE, sep = "\t", row.names = 1,
                               check.names = FALSE, quote = '"', comment.char = ""))
  list(pheno = pdata, expr = expr)
}

reported_sex <- function(pdata) {
  cols <- grep("^characteristics_ch1", names(pdata), value = TRUE)
  hit  <- cols[vapply(cols, function(c) any(grepl("^sex:", pdata[[c]], ignore.case = TRUE)), TRUE)]
  if (!length(hit)) return(rep(NA_character_, nrow(pdata)))
  v <- tolower(trimws(sub("^sex:\\s*", "", pdata[[hit[1]]], ignore.case = TRUE)))
  ifelse(v %in% c("m", "male"), "M", ifelse(v %in% c("f", "female"), "F", NA_character_))
}

sex_probes <- function(idtype, probes) {
  if (idtype == "symbol") {
    hit <- intersect(names(ENTREZ_SEX), probes)
    return(setNames(hit, hit))
  }
  if (idtype == "entrez") {   # "<entrez>_at", or bare "<entrez>" (GSE35488)
    want <- c(paste0(ENTREZ_SEX, "_at"), ENTREZ_SEX)
    hit  <- intersect(want, probes)
    return(setNames(rep(names(ENTREZ_SEX), 2)[match(hit, want)], hit))
  }
  db  <- getExportedValue(idtype, idtype)   # ChipDb keeps the .db suffix
  sym <- suppressMessages(AnnotationDbi::mapIds(db, keys = probes, column = "SYMBOL",
                                                keytype = "PROBEID", multiVals = "first"))
  sym <- sym[!is.na(sym) & sym %in% names(ENTREZ_SEX)]
  setNames(unname(sym), names(sym))
}

# Score = mean z(Y genes) - z(XIST), within cohort; 2-means split. Accepted only
# when the clusters are well separated, so a single-sex cohort is not forced
# into two groups.
infer_sex <- function(expr, idtype) {
  map <- sex_probes(idtype, rownames(expr))
  y   <- intersect(cfg$y_genes, map)
  has_x <- cfg$x_gene %in% map
  if (!has_x && !length(y))
    return(list(sex = rep(NA_character_, ncol(expr)), sep = NA_real_,
                note = "sex marker genes not on platform"))
  m <- collapseRows(expr[names(map), , drop = FALSE], rowGroup = unname(map),
                    rowID = names(map), method = "MaxMean")$datETcollapsed
  z <- t(scale(t(m)))
  # Use whichever side is on the platform (e.g. WG-DASL has the Y genes but
  # no XIST); both sides when available.
  score <- (if (length(y)) colMeans(z[y, , drop = FALSE]) else 0) -
           (if (has_x) z[cfg$x_gene, ] else 0)
  km  <- kmeans(score, centers = range(score))
  sep <- abs(diff(km$centers[, 1])) /
         sqrt(mean(tapply(score, km$cluster, function(s) if (length(s) > 1) var(s) else 0)))
  male <- which.max(km$centers[, 1])
  sex  <- ifelse(km$cluster == male, "M", "F")
  if (!is.finite(sep) || sep < cfg$min_sex_sep)
    return(list(sex = rep(NA_character_, ncol(expr)), sep = sep,
                note = "no clear XIST/Y bimodality"))
  list(sex = setNames(sex, colnames(expr)), sep = sep,
       note = paste(c(if (length(y)) sprintf("%d Y genes", length(y)),
                      if (has_x) "XIST"), collapse = " + "))
}

sex_tab <- list(); sex_log <- list()
for (sid in unique(samples$StudyID)) {
  f <- list.files(file.path(cfg$raw_dir, sub("-GPL[0-9]+$", "", sid)),
                  paste0("^", sid, "_series_matrix\\.txt$"), full.names = TRUE)
  sm <- read_sm(f[1])
  if (sid %in% names(EXPR_FILE)) {
    x <- fread(file.path(cfg$raw_dir, sub("-GPL[0-9]+$", "", sid), EXPR_FILE[[sid]]))
    sm$expr <- as.matrix(x[, -1]); rownames(sm$expr) <- x[[1]]
  }
  ex <- if (sid %in% LINEAR) log2(sm$expr + 1) else sm$expr
  ex <- ex[, colnames(ex) %in% samples[StudyID == sid, GSM], drop = FALSE]
  rep_sex <- reported_sex(sm$pheno)
  if (sid == "GSE93798")   # sex is in the title only: "... H12_M_67_1_Con_L"
    rep_sex <- sub("^.* [A-Z]+[0-9]+_([MF])_.*$", "\\1", sm$pheno$title)
  rep_sex <- setNames(rep_sex, sm$pheno$geo_accession)[colnames(ex)]
  stopifnot(all(rep_sex %in% c("M", "F", NA)))
  inf     <- infer_sex(ex, IDTYPE[[sid]])
  inf_sex <- inf$sex[colnames(ex)]
  sex_tab[[sid]] <- data.table(StudyID = sid, GSM = colnames(ex),
                               Reported = unname(rep_sex), Inferred = unname(inf_sex))
  sex_log[[sid]] <- data.table(StudyID = sid, n = ncol(ex), separation = round(inf$sep, 1),
                               inference = inf$note,
                               reported = sum(!is.na(rep_sex)))
}
sex_tab <- rbindlist(sex_tab)
sex_tab[, Sex := fifelse(!is.na(Reported), Reported, Inferred)]
sex_tab[, SexSource := fifelse(!is.na(Reported), "reported",
                       fifelse(!is.na(Inferred), "inferred", "unknown"))]

conc <- sex_tab[!is.na(Reported) & !is.na(Inferred),
                .(n = .N, agree = sum(Reported == Inferred),
                  concordance = round(mean(Reported == Inferred), 3)), by = StudyID]
sex_log <- merge(rbindlist(sex_log), conc[, .(StudyID, agree, concordance)],
                 by = "StudyID", all.x = TRUE)
sex_log <- merge(sex_log, sex_tab[, .(male = sum(Sex == "M", na.rm = TRUE),
                                      female = sum(Sex == "F", na.rm = TRUE),
                                      unknown = sum(is.na(Sex))), by = StudyID],
                 by = "StudyID")
msg("sex per cohort (inference validated against reported sex where both exist):")
print(sex_log)
write_xlsx(sex_log, file.path(cfg$out_dir, "sex_by_cohort.xlsx"))
write_xlsx(sex_tab, file.path(cfg$out_dir, "sex_by_sample.xlsx"))

# ---------------------------------------------------------------------------
# 3. One multi-factor model per gene
# ---------------------------------------------------------------------------
# A GSM shared by GSE30122 and GSE30528/GSE30529 has one sex; take it per
# (StudyID, GSM) so each observation keeps its own row.
long <- merge(long, sex_tab[, .(StudyID, GSM, Sex, SexSource)], by = c("StudyID", "GSM"))
n_unknown <- uniqueN(long[is.na(Sex), .(StudyID, GSM)])
if (n_unknown) msg("excluded from the model, sex unknown: ", n_unknown, " observations")
dat <- long[!is.na(Sex)]
# The model aggregates across datasets: leave out the subsets of GSE30122.
mdat <- dat[!StudyID %in% cfg$subsets]
stopifnot(!anyDuplicated(mdat[Gene == levels(Gene)[1], GSM]))
msg("model: ", uniqueN(mdat$GSM), " samples, ", uniqueN(mdat$Batch), " batch units")

comps <- sort(unique(mdat$Compartment))
build_design <- function(d) {
  X <- model.matrix(~ 0 + factor(Batch), data = d)
  colnames(X) <- sub("^factor\\(Batch\\)", "Batch:", colnames(X))
  dis <- d$Condition == "Disease"
  for (cp in comps) X <- cbind(X, setNames(data.frame(as.numeric(dis & d$Compartment == cp)),
                                           paste0("Disease:", cp)) |> as.matrix())
  if (uniqueN(d$Sex) > 1) X <- cbind(X, `Sex:M` = as.numeric(d$Sex == "M"))
  X
}

contrast_rows <- function(est, V, df) {
  # Disease effect per compartment, the between-compartment differences, sex.
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
  e  <- drop(L %*% est); se <- sqrt(diag(L %*% V %*% t(L)))
  t  <- e / se
  data.table(Term = rownames(L), log2FC = e, SE = se, t = t, df = df,
             p = 2 * pt(-abs(t), df))
}

fit_gene <- function(d) {
  X <- build_design(d)
  stopifnot("design is rank-deficient" = qr(X)$rank == ncol(X))
  y <- d$Expression
  # Step 1: OLS -> residual variance per batch unit.
  r0 <- lm.fit(X, y)$residuals
  v  <- tapply(r0^2, d$Batch, function(r) sum(r) / max(length(r) - 1, 1))
  w  <- 1 / v[d$Batch]
  # Step 2: WLS with those weights.
  fit <- lm(y ~ 0 + X, weights = w)
  est <- setNames(coef(fit), colnames(X))
  V   <- vcov(fit); dimnames(V) <- list(colnames(X), colnames(X))
  contrast_rows(est, V, fit$df.residual)
}

res <- rbindlist(lapply(levels(mdat$Gene), function(g)
  cbind(Gene = g, fit_gene(mdat[Gene == g]))))
res[, Family := sub(":.*", "", Term)]
# BH within each family across genes/compartments: the 15 disease tests are the
# primary family; between-compartment differences and sex are separate.
res[, FDR := p.adjust(p, "BH"), by = Family]
res[, Term := sub("^Disease:", "Disease vs Control | ", Term)]
res[, Term := sub("^Diff:", "Disease effect difference | ", Term)]
res[, Term := sub("^Sex:M", "Sex (M vs F)", Term)]
res <- res[, .(Gene, Family, Term, log2FC = signif(log2FC, 4), SE = signif(SE, 4),
               t = round(t, 3), df, p = signif(p, 4), FDR = signif(FDR, 4))]
msg("multi-factor model results:")
print(res[Family != "Diff"])
write_xlsx(res, file.path(cfg$out_dir, "stage4_multifactor_results.xlsx"))

# ---------------------------------------------------------------------------
# 4. Per-study Disease effects (sex-adjusted where sex varies) — the per-study
#    estimates the pooling strategy (plan §5) will pool. One dataset at a time,
#    so every dataset.
# ---------------------------------------------------------------------------
per_study <- dat[, {
  ok <- uniqueN(Condition) == 2 && min(table(Condition)) >= 2
  if (!ok) NULL else {
    adj <- uniqueN(Sex) > 1 && all(table(Sex) >= 2)
    f <- if (adj) lm(Expression ~ Condition + Sex) else lm(Expression ~ Condition)
    s <- summary(f)$coefficients["ConditionDisease", ]
    .(n_control = sum(Condition == "Control"), n_disease = sum(Condition == "Disease"),
      sex_adjusted = adj, log2FC = signif(s[1], 4), SE = signif(s[2], 4),
      p = signif(s[4], 4))
  }
}, by = .(Gene, Batch, Compartment)]
write_xlsx(per_study, file.path(cfg$out_dir, "stage4_per_study_effects.xlsx"))
msg("per-study effects: ", nrow(per_study), " gene x batch-unit estimates")

# ---------------------------------------------------------------------------
# 5. Figure — pooled Disease effect per gene and compartment, 95% CI
# ---------------------------------------------------------------------------
pd <- res[Family == "Disease"]
pd[, Compartment := sub("Disease vs Control \\| ", "", Term)]
pd[, `:=`(lo = log2FC - qt(0.975, df) * SE, hi = log2FC + qt(0.975, df) * SE)]
pd[, Gene := factor(Gene, levels = rev(cfg$genes))]
pd[, sig := fifelse(FDR < 0.05, "FDR < 0.05", "n.s.")]
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
out <- file.path(cfg$fig_dir, "stage4_disease_effects.png")
ggsave(out, p, width = cfg$width, height = cfg$height, dpi = cfg$dpi)
msg("wrote ", out)

writeLines(capture.output(sessionInfo()), file.path(cfg$out_dir, "stage4_sessionInfo.txt"))
msg("done")
