# Review of PDF Requirements Against the R Implementation

Reviewed: 2026-10-17. Scope: the non-ML requirements in `CKD_Biomarker_Pooled_Analysis_Plan_Final.pdf` (July 2026, version 2.0), the active R scripts, curated-input metadata, and generated results. Section 6 and the ML ROC plots are excluded. The local decisions in `CLAUDE.md` and `AGENTS.md` are treated as documented scope changes, not as evidence that the original PDF was fully satisfied.

## Assessment

The implementation runs reproducibly and correctly calculates the specified Hedges' g, mean differences, DerSimonian–Laird pooling, Q, and I². Gene mapping and the permitted limma batch-correction route are present. However, it is not a complete implementation of the original non-ML protocol, and the inferential results require attention to biological replication and dependence. Successful execution and agreement between methods do not resolve those design issues.

## Findings, in priority order

### 1. High: repeated biological samples are treated as independent, including repetitions beyond the documented 44 GSM overlaps

**Locations:** `Scripts/stage23_panel_boxplots.R:328–337`; `Scripts/pooling_strategy.R:138–153,177–194`; `Scripts/stage4_MF_stat_mod.R:265–312`; `Scripts/mega_analysis.R:137–145`.

The Stage 2/3 hand-off retains GSM and study identifiers but discards the GEO individual identifiers. Downstream summaries use the number of expression rows as the sample size, and the individual-level models have no patient identifier or correlation structure.

The metadata in `Data/Raw_filtered/GSE30122/GSE30122_series_matrix.txt` demonstrates repeated controls within GSE30122 itself:

- Thirteen glomerular control individuals each have two arrays under different GSM accessions, with matching explicit `individual:` labels. For example, individual 33 appears as GSM757001 and GSM758485; individual 76 appears as GSM757003 and GSM758487.
- Twelve tubular control donor numbers also occur twice in the sample titles. For example, donor 17 appears as GSM757024 and GSM758498, and donor 76 as GSM757029 and GSM758503. Some corresponding `individual:` fields are missing, so title-based linkage is needed.
- Some donors also supply both compartments. Their cross-compartment dependence matters for overall pooling and compartment comparisons.

These observations are not necessarily numerically identical: repeated arrays can have different processed values. Distinct GSM accessions therefore do not establish independent biological replication. The reported 850 distinct GSMs cannot be interpreted as 850 independent patients.

Separately, the known GSE30528/GSE30529 overlap contributes 44 GSMs a second time. I confirmed exact equality of their five-gene values with GSE30122. All 13 series are intentionally retained by repository policy; retaining their inventory entries does not remove the need to account for dependence in statistical inference.

**Consequence:** sample sizes, sampling variances, study weights, confidence intervals, heterogeneity tests, and p-values can be distorted. The existing mega-analysis sensitivity excluding the two subset series does not resolve repetition within GSE30122, and the primary meta-analysis has no equivalent overlap sensitivity. Omitting one analysis unit at a time also leaves related observations in the remaining units.

**Action:** build a biological-sample/individual crosswalk while retaining all 13 series in the inventory. Decide how repeated measurements should contribute to estimation, and account for within-person and overlapping-study dependence. Report a dependence-aware primary or sensitivity analysis. Metafor supports sampling covariance matrices and multilevel structures; a study random intercept alone does not identify the same individual appearing in different study labels. See [metafor's multivariate/multilevel documentation](https://wviechtb.github.io/metafor/reference/rma.mv.html).

A temporary, limited sensitivity excluding GSE30528 and GSE30529 illustrates the impact; it does **not** correct the remaining within-series repetition:

| Gene | Current overall g | Without the two subset series | Current unadjusted p | Sensitivity unadjusted p |
|---|---:|---:|---:|---:|
| C4BPB | -0.0756 | -0.0061 | 0.733 | 0.977 |
| COL13A1 | 0.5743 | 0.5857 | 0.000306 | 0.000772 |
| HABP2 | 0.3757 | 0.2490 | 0.0278 | 0.1406 |
| CFH | 0.4136 | 0.3926 | 0.0753 | 0.0969 |
| CFI | 0.2852 | 0.1983 | 0.1358 | 0.3445 |

HABP2 already has BH-FDR 0.0696 in the current overall analysis; the table must not be described as loss of an originally FDR-significant finding.

### 2. High: reference-RNA replicates are counted as eight independent controls

**Locations:** `Scripts/stage23_panel_boxplots.R:179–181`; `Scripts/pooling_strategy.R:138–149,426`; `Scripts/mega_analysis.R:129–134`.

GSE66494 metadata labels its controls as “Normal kidney total RNA.” Descriptions explicitly identify five discovery-set control replicates and three validation-set control replicates. The code acknowledges reference RNA in its control-type sensitivity, but the main calculations still use `n_Control = 8` as if these were eight independent control subjects.

The metadata does not establish eight independent donors. If the arrays repeat the same reference material, their between-array SD is technical variation rather than between-person variation. It is then inappropriate to interpret the calculated SMD and its sampling variance as a standard patient-versus-control effect. The same concern applies to the single-study mega-analysis and the Stage 4 patient-covariate interpretation.

**Action:** verify donor/material provenance and distinguish technical replication from independent biological controls. Until resolved, label this as a biopsy-versus-reference-RNA comparison and report a sensitivity excluding this cohort from biological-control pooling. This is a provenance/design problem, not an `escalc` arithmetic bug.

### 3. Moderate: the original normalization and clinical requirements remain unmet

**Locations:** `Scripts/stage23_panel_boxplots.R:128,192–209,300–313`; `Scripts/stage4_MF_stat_mod.R:174–244,270–278`.

PDF §3 requires raw processing and uniform quantification/normalization. The active pipeline instead reads submitter-processed series matrices and applies a log transform to the Agilent cohort. It does not perform uniform RMA. `stage1_quantification_old.R` contains archived processing code but is retired and is not the source of active inputs. This is an explicit repository scope change; it should remain an acknowledged protocol deviation. I did not run or restore Stage 1.

Stage 4 adjusts for batch, compartment-specific disease effects, and reported or expression-inferred sex. The curated series-matrix characteristics contain no numeric age, eGFR, proteinuria, or %TIF fields. Consequently, the fuller clinical adjustment and PDF §7 Plot 5 (clinical correlations/scatter plots) are unavailable. Plot 5 is not itself an ML analysis, so excluding ML does not by itself exclude this deliverable.

**Action:** describe the results as sex-adjusted, with inferred sex distinguished from reported sex. Obtain clinical data if the original clinical requirements remain in scope; otherwise formally document them as blocked by unavailable inputs. Do not present their omission as a completed requirement.

### 4. Moderate: failures can silently change the required panel or cohort set

**Location:** `Scripts/stage23_panel_boxplots.R:319–326,347,360–364`.

A missing matrix or probe-extraction error produces a warning and skips the cohort. Only two surviving cohorts are required to continue. If a gene is absent from one cohort, intersecting the gene sets removes it from every cohort, again with only a warning. This can produce downstream analysis files that no longer represent all 13 cohorts and all five required genes.

The current run includes every cohort and gene, so this is a reliability defect rather than an explanation for the current estimates. The final summary count assertion catches some missing-cohort cases late, after earlier outputs have already been written, and does not enforce all five genes.

**Action:** require the expected cohort identifiers and complete five-gene coverage before writing hand-off outputs, or implement an explicit, reviewed missingness policy and exclusion manifest.

### 5. Lower: some reporting overstates what the analysis demonstrates

**Locations:** `Scripts/summary_figures.R:34–41,137`; `Scripts/pooling_strategy.R:347–406`; `Scripts/plots_per_dataset.R:258–264,391–395`.

- The method-agreement title calls the routes independent. They reuse the same samples and expression measurements, and both estimate disease contrasts. Their correlation is descriptive method concordance, not independent biological validation.
- PDF Plot 4 asks for indication-stratified expression including healthy controls and AKI. The integrated heatmap instead shows within-study, control-referenced Hedges' g; controls are a reference rather than a displayed expression stratum, and AKI is absent. This is a defensible adaptation to incompatible expression scales, but it is not literal completion of that plot specification. Per-dataset expression heatmaps provide additional descriptive views.
- Per-dataset Wilcoxon annotations use unadjusted p-values across genes and indication contrasts. This does not violate an explicit PDF correction requirement, but those panels should be labeled exploratory or have a defined multiple-testing correction.

## Requirement-by-requirement status

| PDF requirement | Assessment |
|---|---|
| §2: all five target genes | Implemented; all five are present in all 13 current cohorts. |
| §3 Stage 1: uniform raw processing | Retired; active analysis deviates from the original requirement. |
| §3 Stage 2: mapping and MaxMean collapse | Implemented using platform annotation/custom-CDF Entrez mapping and `WGCNA::collapseRows(method = "MaxMean")`. |
| §3 Stage 3: permitted batch correction | Implemented with `limma::removeBatchEffect`; one full-rank study-by-compartment batch design protects disease effects per compartment. |
| §3 Stage 4: multifactor modeling | Implemented in part: compartment-specific disease contrasts and sex; clinical endpoints unavailable; dependence unresolved. |
| §5: Hedges' g/log2FC and variance | Numerically correct for independent biological groups; replication assumptions need correction. |
| §5: DL random effects, Q and I² | Implemented correctly as specified; independence of supplied effects is the central concern. |
| §5: mega-analysis within platform tiers | Implemented; random study/batch intercept and disease slope with heterogeneous residual variances. Illumina and Agilent are appropriately reported as single-study tiers. Dependence limitations remain. |
| §7 Plot 1: PCA/UMAP before/after | Implemented on the five-gene panel; this demonstrates panel-level behavior, not transcriptome-wide harmonization. |
| §7 Plot 2: five forest plots | Implemented, including individual estimates and pooled diamonds. |
| §7 Plot 4: indication expression heatmaps | Adapted to pooled effects; no AKI cohort or displayed control-expression row in the integrated heatmap. |
| §7 Plot 5: clinical scatter/correlation plots | Not implemented; relevant quantitative outcomes absent from curated metadata. |
| §6 and §7 Plot 3: ML and ROC | Excluded from this review as requested. |

## Correct implementation choices worth retaining

The Stage 4 and pooling scripts use **uncorrected** expression, rather than treating batch-corrected observations as raw modeling inputs. Stage 4 includes batch directly, consistent with [limma's guidance](https://rdrr.io/bioc/limma/man/removeBatchEffect.html).

The batch-correction design intentionally avoids aliased compartment main effects and has a rank check. It protects a separate disease contrast for each compartment. This is a reasonable adaptation of the PDF's simpler study-plus-condition specification. It removes compartment level differences, so corrected figures cannot establish absolute glomerular-versus-tubular expression differences. `removeBatchEffect` is explicitly permitted by the PDF; absence of ComBat or empirical-Bayes shrinkage is not independently a compliance failure.

The pooling code uses SD rather than SE as the `escalc` input, applies the Hedges correction through `measure = "SMD"`, uses `method = "DL"`, and reports Q, I², tau², confidence and prediction intervals. These are appropriate implementations of the specified calculations. See the [metafor effect-size documentation](https://wviechtb.github.io/metafor/reference/escalc.html).

## Validation performed and limits

- Parsed all eight R scripts successfully.
- Ran the five active core scripts in dependency order with only their results directory redirected to `/tmp/renorisk_review/Results`.
- All five completed successfully. All **27** generated Excel workbooks matched the existing corresponding workbooks when read back as tables.
- Confirmed 13 cohorts, 14 batch units, five genes, 894 study/GSM observations, and 850 distinct GSM accessions.
- Independently recomputed mean differences, their sampling variances, and exact small-sample-corrected Hedges' g. Maximum absolute discrepancies were below 7e-16.
- Inspected relevant PDF pages and all curated sample-characteristic field names; checked sample overlaps and repeated individual labels/titles directly against local metadata.
- The cohort-exclusion calculations above are exploratory diagnostics, not replacement results. Other unlinked overlaps across cohorts remain possible. Distinct GSMs are not proof of distinct patients.
- `plots_per_dataset.R` was inspected statically; its full figure-generation run was not repeated. Retired Stage 1 and the network connectivity diagnostic were not run. No analysis scripts or existing results were modified.

The first priorities are a biological-replicate crosswalk, dependence-aware inference, and resolution of GSE66494 control provenance. Protocol deviations and plotting gaps should then be recorded explicitly before describing the work as complete non-ML validation.
