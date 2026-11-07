# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

RenoRisk biomarker validation: a PRISMA-style pooled/meta-analysis of the 5-gene panel
**C4BPA, COL13A1, HABP2, CFH, CFI** across public GEO kidney-biopsy transcriptomics
(CKD, DKD, IgAN, FSGS, AKI, healthy controls). The protocol — roles, cohort inventory,
per-stage specification — is in `CKD_Biomarker_Pooled_Analysis_Plan_Final.pdf`;
the stage numbering below comes from it.

R-first (Bioconductor). R 4.5.2 at `/usr/local/bin/Rscript`. Not a git repository.

## Scope — read this first

**Work only from `Data/Raw_filtered/`. `Data/Raw/` is out of scope.**

Stage 1 (raw CEL/FASTQ re-quantification) has been **retired** — it is preserved as
`Scripts/stage1_quantification_old.R` and reads `Data/Raw`. Do not revive it, propose
Stage 1 runs, or reach for the salmon/rsubread/ncbi engines without asking first.

**RNA-Seq is in the same bag as the microarrays (by instruction, 2026-11-03).**
`Data/Raw_filtered/` holds three RNA-Seq series: GSE137570 (GPL23227 only; the older
GPL11154 part was deleted by the user), GSE142025, GSE219185. Their series matrices are
metadata-only; the values are in `<id>_expression.txt` (submitter-normalised log2, gene ×
GSM — not counts, so DESeq2/VST/ComBat-seq do not apply to them; GSE219185 is voom log2-CPM with
Ensembl row IDs). Treatment (user decision, 2026-11-03):
- **GSE142025** (27 DN vs 9 tumor-nephrectomy controls, whole biopsy) is a regular
  Stage 2/3 cohort (`idtype = "symbol"`, `expr_file`) and so enters EVERY analysis:
  batch correction, Stage 4 model, pooling, mega-analysis (its own `RNA-Seq` tier),
  per-dataset plots, summaries, ML inputs. C4BPA is not in its deposit.
- **GSE137570** (24 CKD, no controls; has age/sex/eGFR/%TIF) and **GSE219185** (220
  NEPTUNE tubulointerstitium, FSGS vs MCD, no controls) enter the **ML inputs only**, via
  `Scripts/by-claude/rnaseq_panel.R`. Without controls their batch offset cannot be separated from
  disease and they give no within-study effect.
- GSE219185 shares patients with NEPTUNE GSE108109 (glomerular) — keep them in one LODO
  group. C4BPA is NA for GSE142025 and GSE137570 in the ML tables.
- **GSE141295** (glomeruli, 14 IgAN vs 10 normal cortex of renal-cell-carcinoma kidneys
  = Tumor nephrectomy) and **GSE175759** (tubulointerstitium, 65 disease — IgAN 43, MCD 9,
  MN 4, DN 3, FSGS 3, LN 3 — vs 22 "nephrectomy controls"; per-sample eGFR, unused) are
  regular cohorts in EVERY analysis (by instruction, 2026-11-03). GEO holds their values
  only as raw counts: `<id>_RAW.tar` was downloaded into the folder and extracted to
  `<id>_RAW/`; **`Scripts/by-claude/rnaseq_counts_vst.R` (run first) writes `<id>_expression.txt`**
  — Ensembl → symbol (counts summed per symbol), then **DESeq2 VST** (by instruction,
  2026-11-04; replaced TMM log2-CPM): design `~ condition` (the cohort's label field),
  genes with `rowSums >= 10`, `vst(blind = FALSE)` — the same steps as the per-dataset
  RNA-Seq boxplot branch (`Examples/boxplot_rr_dnaseq.r`), after ComBat-seq on the two
  series together (see Invariants). It also writes `<id>_counts.txt`
  (symbol-summed counts); registry flag `counts = TRUE` sends their per-dataset boxplots
  to the DESeq2 branch (Wald padj), which plots the `<id>_expression.txt` values — VST
  shifts all genes by a near-constant offset (0.03 / 0.16) when only the factor-level
  names change, so it is not recomputed for the plots. Normalisation only, not re-quantification. GSE175759's 3
  submitter-flagged technical outliers (all IgAN) are left out of the file; loaders drop
  series-matrix samples that have no values. All 5 panel genes pass the count filter in both. Only these two series
  have raw counts: GSE142025 (normalised log2), GSE137570 (quantile-normalised) and
  GSE219185 (voom) cannot take DESeq2. QC: `Results/RNASeq/counts_vst_qc.xlsx`.

**Eight series were added to `Data/Raw_filtered/` on 2026-11-03; three are included (by
instruction, same day)** — assessment in `Raw_filtered_inclusion_review.docx`:
- Included in every analysis (regular Stage 2/3 cohorts): **GSE93798** (glomeruli, 20 IgAN
  vs 22 controls, custom CDF v19 as GSE104948; sex only in the sample title, read as
  reported — `title_sex`; controls "Con_L" = living donors, confirmed by the user 2026-11-04), **GSE35488** (tubulointerstitium,
  25 IgAN vs 6 living donors, Entrez CDF v10 — **bare Entrez IDs, no `_at`**; the `entrez`
  route accepts both), **GSE111154** (post-mortem kidney cortex, 4 early DN vs 4
  non-diabetic donors, HTA 2.0; compartment "Whole biopsy", control type "Deceased donor").
  Their `.txt.gz` series matrices were decompressed beside the originals (`gunzip -k`).
- Duplicate folders are renamed `<GSE>-dup` (by instruction, 2026-11-03: marked for removal):
  `GSE35487-dup`, `GSE37460-dup`, `GSE30528-dup`, `GSE30529-dup`. No script reads them.
- Folders excluded for another reason are renamed `<GSE>-excl` (by instruction, 2026-11-04):
  `GSE1009-excl` — glomeruli, 2 normal vs 2 DN kidneys, each RNA hybridised in duplicate
  (Baelde 2004, PMID 15042541; GEO holds 6 of the 8 arrays: Control/Diabetes 1a+1b = one
  patient each). 2 patients per arm is below the pipeline minimum of 3 — the only reason;
  platform (HG-U95Av2, `hgu95av2.db` not installed) and linear MAS5 values were fixable.
- Excluded (PRISMA reasons in `cfg$excluded`): GSE35487 (same 31 arrays as GSE35488,
  other CDF), GSE37460 (GPL11670 = the 18 GSE104948 living donors, same codes; GPL14663 =
  40/51 already in GSE99339 under the same `Hx-Glom-...` codes), GSE1009 (HG-U95Av2, 2
  patients per group hybridised twice, MAS5 linear). GSE141295/GSE175759 were first
  excluded (counts only) and then included — see above.
- GSE137570 also gained `_clinical.txt` (same fields as its series matrix) and a
  normalised-counts xlsx (not read).

## Pipeline stages

| Stage | Status | Where |
|---|---|---|
| 1 — raw processing & uniform quantification | **retired** | `stage1_quantification_old.R` |
| 2 — target gene mapping | implemented | `stage23_panel_boxplots.R` |
| 3 — batch correction | implemented | `stage23_panel_boxplots.R` |
| 4 — multi-factor statistical modelling | implemented | `stage4_MF_stat_mod.R` |
| RNA-Seq raw counts → DESeq2 VST input (GSE141295, GSE175759) | implemented | `by-claude/rnaseq_counts_vst.R` → `Data/Raw_filtered/<id>/<id>_expression.txt` (run first; not in `all.R`) |
| RNA-Seq series without controls → ML inputs | implemented | `by-claude/rnaseq_panel.R` → `Results/RNASeq/` |
| pooling strategy (plan §5; not a stage) | implemented | `pooling_strategy.R` |
| mega-analysis within platform tiers (plan §5) | implemented | `mega_analysis.R` |
| summary figures (method agreement, PRISMA flow) | implemented | `by-claude/summary_figures.R` |
| ML input workbook | implemented | `create_ML_train_data.R` (+ `by-claude/rnaseq_panel.R`); also Part I of `all.R` (same output) |
| ML input, all common genes | implemented | `all2.R` → `Data/ML_input/all2.xlsx` (1322 × 8,314 genes incl. all 5 RNA-Seq series; C4BPA is NOT among them — absent from GSE142025/GSE137570; self-contained, reads the series matrices + RNA-Seq `_expression.txt` — keep its registry/rules in sync with `all.R`) |
| ML input, ComBat-corrected | implemented | `create_ML_corr_data.R` → `Data/ML_input/all_corr.xlsx` (run after `create_ML_train_data.R`) |
| LODO ML (plan §6) | not implemented | — |

Stage 4 as the plan defines it is only the multi-factor DE model. Still to do, from
plan §6: LODO (`LeaveOneGroupOut`, grouped by study) logistic/RF/SVM classification and
ElasticNet regression on eGFR and %TIF (Python/scikit-learn, `random_state=42`) — note
no microarray cohort carries eGFR or %TIF (see Stage 4 below); only the RNA-Seq series
GSE137570 (GPL23227, n = 24; eGFR and %TIF) and GSE175759 (n = 87; eGFR) carry them, and
neither is wired in.

**ML dataset inclusion rule (user instruction, 2026-10-19): exclude the GSE30528
and GSE30529 series when constructing training/validation data.** Their samples
remain represented through GSE30122 only; do not remove those GSMs from the parent
series. Since 2026-11-03 both are dropped from every analysis (folders `-dup`), so the
rule holds trivially; `create_ML_train_data.R` keeps the filter as a guard.

`Scripts/create_ML_train_data.R` writes `Data/ML_input/all.xlsx` (by instruction, 2026-10-30; readers:
`train_ML_RR.ipynb`, `DNN-sample/redemption_train.py`, `visualize_ML_data.py`) from the
Stage23 uncorrected panel and indication tables plus Stage4 `sex_by_sample.xlsx`, plus
`Results/RNASeq/` (GSE137570, GSE219185): 1322 GSM rows (931 microarray + 391 RNA-Seq;
1077 Disease / 245 Control, 19 studies), five expression features (C4BPA NA for GSE142025
and GSE137570), metadata, `Platform` (`Microarray` / `RNA-Seq`), `Sex` (`M` / `F`; missing
values remain NA; by instruction, 2026-10-22 — no `SexSource` column), and
`Disease_Label` (Control = 0, Disease = 1). Reported sex takes priority, otherwise
inferred; per-sample provenance stays in Stage4 `sex_by_sample.xlsx`. Sex is an
optional categorical predictor. No additional
scaling, cross-study batch correction, overlap-based exclusions or train/test
splits are applied. A README sheet records provenance and overlap limitations;
`ML_training_sessionInfo.txt` sits beside the workbook. The script verifies the
written values by reading them back. Run after Stage23, Stage4 and `by-claude/rnaseq_panel.R`.
It also writes `all_bn.xlsx` (by instruction, 2026-10-20): same rows, columns and
metadata, expression = the Stage23 `removeBatchEffect` model refitted gene by gene over
all 1322 rows (so the disease-only RNA-Seq series are corrected too); values differ from
`panel_harmonised_long.xlsx` (fitted on 1078 rows; with the plan's `~ StudyID + Condition`
design the difference is no longer one constant per gene). Its README flags
label leakage — the correction design protected Disease, so CV on it is optimistic.

`Scripts/create_ML_corr_data.R` writes `Data/ML_input/all_corr.xlsx` (by instruction,
2026-11-02): `all.xlsx` with `sva::ComBat` (parametric EB, location + scale) over all 1322
samples; batch = `StudyID` (19 studies), `mod` = intercept + Disease — plan §3
"~ StudyID + Condition" (by instruction, 2026-11-04; it replaced batch = 20 study ×
compartment units with a Disease effect per compartment). `batch_r2` is still measured
against the 20 `Batch` units.
Only the five gene columns change. Sheets: `training_data`, `README`, `batch_r2`
(per-gene variance share explained by batch: ~0.90–0.98 before, ≤ 0.02 after).
ComBat cannot fit a batch with no values for a gene, so it runs in two passes: the
complete genes over all rows, then C4BPA (with all 5 genes, for the EB priors) over the
1262 rows that have it.
Same label-leakage caveat as `all_bn.xlsx`, and with 5 genes the EB priors are weak.
ComBat-seq runs earlier, on the two count-deposited series only (`by-claude/rnaseq_counts_vst.R`); DerSimonian–Laird pooling
yields per-study summaries, not sample values, so it is not an ML-input transform.

### Pooling strategy (plan §5) — `Scripts/pooling_strategy.R`

A random-effects meta-analysis; "pooling strategy" is the plan's §5 name for it (renamed
from `meta_analysis.R` by instruction, 2026-10-09). Plan §3 defines only Stages 1–4.

Reads `panel_uncorrected_long.xlsx` (effects are within-study, so uncorrected values
are correct). Units = the Stage 2/3 batch units: effect sizes for, and pooling over, all
18 (`Pooled` column; 17 for C4BPA, absent from GSE142025). Glomerular k = 10, Whole
biopsy k = 4 (GSE116626, GSE66494, GSE142025, GSE111154), Tubulointerstitial k = 4
(GSE30122 tubuli, GSE104954, GSE35488, GSE175759); GSE142025 also
enters the indication heatmap (DKD) and is typed Tumor nephrectomy in the control-type
table (whole biopsy, so not in the glomerular test). Hedges' g (`escalc` SMD; primary,
scale-free — matters because four cohorts are per-gene centred) and log2FC (MD;
secondary), pooled by `rma(method = "DL")`: overall, per compartment, plus a
`mods = ~ Compartment` test; Q, I², tau², prediction intervals; `leave1out`. One
metafor warning "SMD values unusually large" is expected: CFH in GSE30122 tubuli is a
real g = 2.3, not an SE/SD mix-up.

### Stage 4 — `Scripts/stage4_MF_stat_mod.R`

Reads `Results/Stage23/panel_uncorrected_long.xlsx` (Stage 2/3 values **before**
`removeBatchEffect` — limma advises against modelling corrected values). Per gene, one
WLS model: `Expression ~ 0 + Batch + Disease:Compartment + Sex`, weights = 1/residual
variance per batch unit (spread differs up to ~6.7× between studies). Compartment main
effects are aliased with batch, so it reports the Disease effect per compartment and the
between-compartment differences of that effect. BH within family (Disease / Diff / Sex).

**No eGFR, %TIF, proteinuria or age exists in any of the 16 microarray series matrices**
(the RNA-Seq GSE137570 has age, sex, eGFR, %TIF — unused) — the plan
annex says otherwise for ERCB and GSE66494, but those values live in Nephroseq, not GEO.
Sex is the only patient covariate: reported for GSE116626, GSE129973, GSE96804; inferred
elsewhere from XIST vs Y genes (validated 40/40 and 60/61; rejected for GSE116626, whose
WG-DASL chip lacks XIST, but it reports sex). Stage 4 must run after Stage 2/3.

## Current pipeline — `Scripts/stage23_panel_boxplots.R`

Self-contained Rscript: `cfg` block, a `COHORTS` registry, helpers, driver. Runs in ~30s.
Reads the series matrices directly (metadata + expression table); no GPL fetch.

Flow: parse → map probes to the panel → MaxMean collapse (`WGCNA::collapseRows`) → log2 the one linear cohort →
**one** `limma::removeBatchEffect` model over all samples (plan §3 "~ StudyID + Condition":
batch = StudyID, design = intercept + Disease) → plot.

To change genes, edit `cfg$genes` **and** `ENTREZ` (the custom-CDF cohorts address probes
as `<entrez>_at`, so symbol alone is not enough). Everything else is driven off those.

### Invariants — each encodes a correctness decision

- **Duplicate datasets are dropped; their folders are renamed `<GSE>-dup`** (by
  instruction, 2026-11-03; it replaced the 2026-10-19 rule "GSE30528/GSE30529 in
  single-dataset calculations only"). GSE30122 fully contains GSE30528 (22/22 GSMs) and
  GSE30529 (22/22); both are now out of every analysis, as are GSE35487 and GSE37460.
  No contained dataset is counted: 17 series (14 microarray + 3 RNA-Seq), 18 units, 1078
  samples in every analysis. **Patient-level overlap without shared arrays (the same
  donor biopsied twice, deposited twice) cannot be detected from these files** — accepted
  by the user (2026-11-04) as a standing limitation. The "fully contained → single-dataset use only" guard is kept but
  empty: `cfg$subsets <- character(0)` in each script (`all.R`: registry field
  `contained_in` -> `SUBSETS`). **Partially overlapping datasets stay as they are — do
  not fix, trim or modify them** (by instruction, 2026-11-03): only a dataset whose data
  is wholly contained in another is excluded. Known partial overlaps, kept: GSE108109 /
  GSE200828 (96 shared reprocessed samples; Results/SampleOverlap/README.md), GSE121233 /
  GSE129973 (candidate). GSE37460 (40/51 arrays in GSE99339, its 18 LDs in GSE104948) was
  excluded by explicit instruction.
- **One batch-correction model over all aggregated samples** (by instruction, 2026-10-08; it
  replaced three per-compartment fits). Compartment is nested inside study for most
  cohorts, so compartment main effects are NOT in the design — they are aliased with
  batch, and `removeBatchEffect` silently zeroes aliased batch coefficients
  (`beta[is.na(beta)] <- 0`), leaving those studies uncorrected. A `qr()` rank check
  guards this. **Model = the plan's `~ StudyID + Condition`** (by instruction, 2026-11-04;
  it replaced batch = study × compartment unit + a Disease effect per compartment):
  batch = StudyID (17), design = intercept + Disease, here, in `all_bn.xlsx` and in
  the `all_corr.xlsx` ComBat. The `Batch` column (study × compartment, GSE30122 split
  on `source_name_ch1`, 18 units) stays the unit of Stage 4, pooling, mega-analysis and
  the PCA/UMAP batch metrics. Consequence: cross-compartment level differences are
  removed except inside GSE30122 (one batch, both compartments) — compare within
  compartment only.
- **ComBat-seq on raw counts** (plan §3; by instruction, 2026-11-04): only GSE141295 and
  GSE175759 have counts. `by-claude/rnaseq_counts_vst.R` runs `sva::ComBat_seq` over both on their
  34,130 common genes (batch = StudyID, group = Disease/Control) before the DESeq2 VST,
  so their values — also the "uncorrected" Stage 4 / pooling input and the per-dataset
  plots — are ComBat-seq adjusted. Stage 3 then corrects them again with the rest.
- **A panel gene missing from a deposit stays in the panel** (C4BPA in GSE142025). That
  cohort is NA for it; `removeBatchEffect` is fitted gene by gene on the samples that have
  the gene (each with its own rank check — a batch unit with no values would otherwise be
  aliased and silently zeroed); PCA/UMAP use the genes on every cohort (4). Refactor
  verified: without GSE142025 every Stage23 table matched the previous run (≤ 1e-12).
- **`removeBatchEffect` fixes location, not spread.** Residual per-gene SD still differs
  up to ~6.7× across studies (`residual_sd_ratio.xlsx`). Measure this per gene — pooling
  genes is dominated by between-gene variance and looks falsely uniform.
- **Multi-probe genes are collapsed with `WGCNA::collapseRows(method = "MaxMean")`**
  (by instruction, 2026-10-08). It replaced an inline `which.max(rowMeans)` collapse;
  all Stage23 outputs were verified byte-identical across the switch.
- Annotation objects are exported under the **full package name including `.db`** —
  `getExportedValue("hgu133a2.db", "hgu133a2.db")`. Stripping the suffix yields a
  function and a confusing `mapIds` dispatch error.

## The curated set — facts that cost a full survey to derive

14 microarray cohorts used from `Data/Raw_filtered/` (11 original + GSE93798, GSE35488,
GSE111154), plus 5 RNA-Seq series described under Scope (GSE142025, GSE141295, GSE175759
full cohorts, ID route = gene symbol; GSE137570, GSE219185 ML only), plus 5 excluded folders
(see Scope; the 4 duplicates carry `-dup`, GSE1009 carries `-excl`). Seven array platforms, four probe-ID
namespaces, three incompatible value scales. Nothing here is recorded in the data itself.

| Cohorts | ID route | Scale as deposited |
|---|---|---|
| GSE104948-GPL22945, GSE104954-GPL22945, GSE99339-GPL19184, GSE108109, GSE200828, GSE93798 | Entrez `<id>_at` (custom CDF) | log2 absolute |
| GSE35488 | bare Entrez `<id>` (custom CDF v10), same `entrez` route | log2 absolute |
| GSE121233, GSE129973, GSE96804, GSE111154 | `hta20transcriptcluster.db` | log2 absolute (GSE111154: Partek, ~4 units higher) |
| GSE30122 | `hgu133a2.db` | **per-gene centred** (row median exactly 0) |
| GSE116626 (Illumina BeadChip) | `illuminaHumanv4.db` | **per-gene centred** |
| GSE66494 (Agilent) | `hgug4112a.db` | **linear — the only cohort needing log2** |

Not all cohorts are Affymetrix: 12 Affymetrix series, 1 Illumina, 1 Agilent. The per-gene-centred
cohorts have lost absolute expression irrecoverably — any figure sharing a Y axis with
them is relative, and must say so.

**Condition labels hide in a different field per cohort**: `diagnosis` with a blank value
meaning living donor (ERCB; its 3+3 "Tumor Nephrectomy" samples are **also Control**, by
instruction 2026-10-09, matching GSE99339/GSE129973/GSE96804 where tumor nephrectomy is
the only control group — ERCB controls are 21 per compartment, not 18), `subject status` (GSE116626), `sample type` (GSE66494),
`description` (GSE99339), `tissue1` (GSE141295), `diagnosis` (GSE175759), `source_name_ch1`
(the rest, incl. GSE93798, GSE35488, GSE111154). The registry's `cond_field` /
`control` columns encode this; `derive_condition()` applies it.

Compartment strata in the aggregated analyses: Glomerular (10 units, 680 samples),
Tubulointerstitial (4, 212), Whole biopsy (4, 186 — incl. GSE142025, GSE111154).

## Running

```bash
Rscript Scripts/by-claude/rnaseq_counts_vst.R  # first: GSE141295/GSE175759 counts -> <id>_expression.txt
Rscript Scripts/stage23_panel_boxplots.R   # Stage 2/3
Rscript Scripts/stage4_MF_stat_mod.R       # Stage 4 (reads Stage 2/3 output)
Rscript Scripts/pooling_strategy.R         # plan §5 (reads Stage 2/3 output)
Rscript Scripts/mega_analysis.R            # plan §5 mega-analysis (reads Stage 2/3 + pooling output)
Rscript Scripts/by-claude/summary_figures.R          # last: reads Stage 2/3, Stage 4 and pooling output
Rscript Scripts/by-claude/rnaseq_panel.R             # GSE137570 / GSE219185 for the ML inputs
Rscript Scripts/create_ML_train_data.R     # all.xlsx + all_bn.xlsx (after Stage 2/3, Stage 4, rnaseq_panel)
Rscript Scripts/create_ML_corr_data.R      # all_corr.xlsx (after create_ML_train_data.R)
```

No arguments. Requires `limma`, `ggplot2`, `data.table`, `AnnotationDbi`, `WGCNA`, `uwot` (UMAP), `writexl`, `readxl`, `metafor` (pooling strategy), `nlme` (mega-analysis), `sva` (ComBat, `create_ML_corr_data.R`), `DESeq2` (`by-claude/rnaseq_counts_vst.R`) plus the four
platform annotation packages (`hgu133a2.db`, `hta20transcriptcluster.db`,
`illuminaHumanv4.db`, `hgug4112a.db`) and `org.Hs.eg.db` (GSE219185 Ensembl IDs) — all installed. Unlike the retired Stage 1, this
script does **not** auto-install and needs no network.

No test suite. Verify a change by re-running and checking: `cohort_qc.xlsx` lists 17
cohorts (5 genes each, GSE142025 4), the log says "one model: corrected 1078 samples across 17 batch units = studies",
and `panel_boxplot_per_cohort.png` shows cohorts aligned at a common level per gene across
all compartments (that panel is the did-the-correction-work check).

## Data layout

- `Data/ML_input/` — the generated ML inputs (the only generated files kept under `Data/`):
  `all.xlsx` (5-gene panel, by `create_ML_train_data.R` / `all.R`) and `all2.xlsx` (all common
  genes, by `all2.R`; read by `train_ML_ALL.ipynb`) and `all_corr.xlsx` (`all.xlsx` after
  ComBat, by `create_ML_corr_data.R`; label leakage — not for unbiased CV). `all_bn.xlsx`
  stays in `Results/ML_training/`.
- `Data/Raw_filtered/<GSE>/` — **the input.** 24 series folders (19 used: 14 microarray +
  5 RNA-Seq; 5 excluded: 4 duplicates renamed `<GSE>-dup`, GSE1009 renamed `-excl`), one platform per used GSE (ERCB → GPL22945, GSE99339 → GPL19184,
  GSE137570 → GPL23227); scripts read `<id>_series_matrix.txt`, not `.gz`; the 5 RNA-Seq folders
  add `<id>_expression.txt` (read; for GSE141295/GSE175759 derived from `<id>_RAW/` counts)
  and supplementary files (not read).
- `Results/` — **everything the pipeline generates**, one subfolder per output group
  (each script sets `cfg$results_dir`). `Data/` holds inputs only. **Each script writes its
  tables, figures and `<name>_sessionInfo.txt` into one folder** (`cfg$fig_dir` =
  `cfg$out_dir`; by instruction, 2026-10-09 — there is no shared figure folder any more).
  - `Results/Stage23/` — by Stage 2/3. Figures: `panel_boxplot_by_compartment.png` (the
    deliverable figure), `panel_boxplot_per_cohort.png` (QC), `pca_before_after_correction.png`,
    `umap_before_after_correction.png` (plan §7 Plot 1). Tables: `panel_harmonised_long.xlsx`,
    `panel_uncorrected_long.xlsx` (Stage 4 and pooling-strategy input),
    `cohort_qc.xlsx`, `residual_sd_by_study_gene.xlsx`, `residual_sd_ratio.xlsx`,
    `sample_indication.xlsx` (specific diagnosis per sample, from the
    `INDICATIONS` patterns; kept out of the 7-column hand-off file), `pca_batch_r2.xlsx` (share of PC1/PC2 explained by batch unit and by Condition,
    before vs after correction), `umap_batch_silhouette.xlsx` (mean silhouette width of
    the batch units on the UMAP — understates the "before" batch effect, because cohorts on
    the same platform cluster together and GSE30122/30528/30529 are the same samples).
  - `Results/Stage4/` — **all** Stage 4 output (tables and figure):
    `stage4_multifactor_results.xlsx`, `stage4_per_study_effects.xlsx`, `sex_by_cohort.xlsx`,
    `sex_by_sample.xlsx`, `stage4_disease_effects.png`.
  - `Results/PoolingStrategy/` — by `pooling_strategy.R`. Figures: `pooling_forest_<GENE>.png`
    (plan §7 Plot 2), `pooling_summary_forest.png` (all genes × compartments, CI + prediction
    interval), `pooling_indication_heatmap.png` (plan §7 Plot 4),
    `pooling_control_type_sensitivity.png`. Tables: `pooling_study_effects.xlsx`,
    `pooling_estimates.xlsx`, `pooling_compartment_moderator.xlsx`, `pooling_leave1out.xlsx`,
    `pooling_by_indication.xlsx` (Hedges' g per indication vs in-study controls, DL-pooled),
    `pooling_control_type_{study_effects,glomerular,test}.xlsx` — sensitivity to the control
    tissue: Living donor (ERCB restricted to its LDs, GSE108109, GSE200828, GSE93798,
    GSE116626, GSE35488) vs Tumor nephrectomy (GSE99339, GSE129973, GSE96804, GSE121233 —
    per its overall design — GSE141295, GSE142025); GSE30122 = "Not stated", GSE66494 =
    reference RNA replicates, GSE111154 = Deceased donor, GSE175759 = Nephrectomy. Tested in
    glomeruli only (k = 4 vs 5): **CFI differs** (TN − LD g = 1.27, FDR 0.035); the other
    genes do not (FDR ≥ 0.37). (Before GSE93798/GSE141295 joined: k = 3 vs 4, all FDR ≥ 0.32.)
  - `Results/MegaAnalysis/` — by `mega_analysis.R` (plan §5, "mega-analysis within platform
    tiers"). One `nlme::lme` per tier × gene on uncorrected log2 values: fixed Disease effect
    per compartment, random study intercept + disease slope (`pdDiag`), residual variance per
    study (`varIdent`). Tiers: Affymetrix (12 series, 13 units: all array glomerular + tubulointerstitial
    units + GSE111154 whole tissue — the only tier that pools), Illumina (GSE116626), Agilent
    (GSE66494) and RNA-Seq (GSE142025 whole biopsy, no C4BPA; GSE141295 glomerular; GSE175759
    tubulointerstitial). **A compartment with one study in its tier gets that study's own
    `lm`** — a joint lme would borrow the disease-slope variance from other compartments. Tables: `mega_estimates.xlsx`, `mega_vs_meta.xlsx`, `mega_vs_meta_agreement.xlsx`, `mega_tiers.xlsx`; figure
    `mega_vs_meta.png`. A z-scored (standardized) scale was tried and dropped: z-scoring by
    each study's total SD collapses the between-study disease variance to 0.
  - `Results/RNASeq/` — by `by-claude/rnaseq_panel.R` (GSE137570, GSE219185; ML inputs only):
    `panel_rnaseq_long.xlsx` (same 7 columns as `panel_uncorrected_long`),
    `rnaseq_samples.xlsx` (+ Indication, Label, Sex), `rnaseq_qc.xlsx`.
  - `Results/Summary/` — by `by-claude/summary_figures.R`. Figures: `method_agreement.png` (Stage 4 vs
    pooling), `prisma_flow.png`. Tables: `method_agreement.xlsx`,
    `method_agreement_stats.xlsx`, `prisma_counts.xlsx`. PRISMA: 16 datasets listed in the
    curated dataset inventory (an Excel sheet kept outside the project; same 16 accessions
    as plan Annex 2) + 8 added to Raw_filtered on 2026-11-03 = 24 → 7 excluded (no control
    arm: GSE137570, GSE219185 — used in the ML inputs only; duplicates: GSE35487, GSE37460,
    GSE30528, GSE30529; too few patients: GSE1009) → 17 included (14 microarray + 3 RNA-Seq).
    Adding 5 series beyond Annex 2 (GSE93798, GSE35488, GSE111154, GSE141295, GSE175759) is a
    protocol deviation — reports must say so. The plan text's "32 datasets" is NOT the selection count — do not use it.
  - `Results/SingleDatasetPlots/<GSE>/` — by `plots_per_dataset.R`: per-gene boxplots
    and two heatmaps per dataset — `_heatmap_top50.png` (top-50 variable probes) and
    `_heatmap_<n>genes.png` (the panel genes found: 5, or 4 for GSE142025); plus
    `boxplots_summary.xlsx` (`missing_genes` column). GSE142025 is plotted like an array
    (groups Control / Early_DN / Advanced_DN), not through the DESeq2 branch. GSE141295 and GSE175759 go
    through the DESeq2 branch (VST values, DESeq2 Wald padj vs the reference arm).
- `Data/Raw/` — out of scope (see above). Never commit; large and re-downloadable.
- `.claude/skills/` — 177 vendored skills from K-Dense-AI/scientific-agent-skills (MIT).
  Project-scoped, Python-leaning; `bulk-rnaseq` and `pydeseq2` do not fit this R-first,
  microarray-only pipeline.

## Conventions

- **`Scripts/by-claude/`** (by instruction, 2026-11-04) holds the scripts Claude wrote without an
  explicit "create this script" request: `rnaseq_counts_vst.R`, `rnaseq_panel.R`,
  `summary_figures.R`, `audit_expression_overlap.R`. Run them as `Rscript Scripts/by-claude/<name>`;
  they use the absolute `cfg$project_dir`, so their inputs and outputs are unchanged. Put any new
  script of that kind there too.

- **`all.R` and the modular scripts change together, always, in both directions** (by
  instruction, 2026-11-03). `all.R` is a single-file mirror: Part B/D =
  `stage23_panel_boxplots.R`, C = `plots_per_dataset.R`, E = `stage4_MF_stat_mod.R`,
  F = `pooling_strategy.R`, G = `mega_analysis.R`, H = `by-claude/summary_figures.R`,
  I = `create_ML_train_data.R` + `by-claude/rnaseq_panel.R` (inline; writes the same `Results/RNASeq/`
  tables). Any edit to one side is made on the other in the same change. Verify by running
  the modules, snapshotting `Results/` + `Data/ML_input/`, running `all.R`, and comparing
  every `.xlsx` it rewrites (`all.equal`, tolerance 1e-10): on 2026-11-03 all 33 matched (again after adding GSE93798/GSE35488/GSE111154).
  Not in `all.R`: `create_ML_corr_data.R` (`all_corr.xlsx`), `all2.R` (`all2.xlsx`) and
  `by-claude/rnaseq_counts_vst.R` (its `<id>_expression.txt` files are inputs `all.R` reads).
- `cfg$project_dir` is a hard-coded absolute path (`~/Striki/Quant/quant.2/MetaAnalysis`);
  new scripts follow the same `cfg` list + derived paths pattern.
- **All tabular data the pipeline writes or reads is `.xlsx`** (by instruction, 2026-10-09):
  `writexl::write_xlsx` to write, `readxl::read_excel` + `as.data.table()` to read, one
  workbook per table. Readers pass explicit `col_types` — never let readxl guess. Only
  the GEO series matrices (`.txt` inputs) are read as text. Missing values are true `NA`
  (empty cell), not the `""` the old TSV round trip produced.
- **Input provenance is a standing caveat.** These are submitter-processed values from
  five different pipelines (GeneSpring, Expression Console, a "modified Affy package",
  Illumina quantile, Agilent FeatureExtraction). Stage 1 existed to bypass exactly this;
  with it retired, every figure carries that limitation and should state it.
- Every stage writes provenance next to its outputs (`*_sessionInfo.txt`, QC tables) —
  the Methods section is assembled from these.
- `set.seed(1)` at the top.
