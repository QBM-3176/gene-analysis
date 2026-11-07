# RenoRisk implementation audit

Date: 2026-10-17. Scope: PDF requirements other than the ML stage, active R scripts, curated series matrices, and existing results. No analysis scripts or existing analysis outputs were modified. This review treats the repository's documented decisions as authorized scope decisions, while separately assessing compliance with the original PDF and statistical validity.

## Verdict

The active pipeline is reproducible and implements the main mapping, batch-correction, differential-expression, random-effects meta-analysis, and platform-tier mega-analysis operations. It is **not a complete implementation of the original non-ML protocol**, and its current inferential results need qualification. The most consequential issues concern dependent observations and reference controls, rather than syntax or incorrect calls to `metafor`.

## Requirements comparison

| PDF requirement | Assessment | Evidence and limitation |
| --- | --- | --- |
| Section 3, Stage 1: uniform raw quantification, Affymetrix RMA | Outside the active implementation | Stage 1 is retired by repository instruction. Active scripts use submitter-processed matrices, some already batch corrected. This is an explicit departure from the PDF, not a missing invocation to restore automatically. |
| Stage 2: all five genes, WGCNA maximum-mean probe collapse | Implemented | `stage23_panel_boxplots.R:282` maps custom-CDF Entrez IDs or platform annotations; line 309 calls `collapseRows(method = "MaxMean")`. All five genes are present in all 13 active series. |
| Stage 3: continuous-expression batch correction with study and condition | Implemented with a documented design adaptation | Line 403 calls `limma::removeBatchEffect`; batch is study × compartment, with protected compartment-specific disease effects and a full-rank check. No RNA-seq exists in the active curated inputs, so ComBat-seq is inapplicable to those inputs. |
| Stage 4: compartment and patient clinical covariates | Partial | Batch terms absorb compartment baselines; disease effects and their differences are estimated per compartment. Sex is included, mostly inferred. eGFR, fibrosis, age, proteinuria, and CKD stage are unavailable in the curated sample characteristics and are not modelled. |
| Section 5: per-study Hedges' g or log2FC and sampling variance | Implemented, with validity issues below | `pooling_strategy.R:138` computes within-unit means, SDs, and counts, then `escalc("SMD")` and `escalc("MD")`. |
| Section 5: DL random effects, Q and I² | Implemented | `pooling_strategy.R:177` uses `rma(..., method = "DL")`; Q, I², tau², confidence intervals, prediction intervals, and leave-one-out results are exported. Standard independence assumptions are violated. |
| Section 5: mega-analysis within platform tiers | Implemented with documented modelling choices | `mega_analysis.R:127` fits sample-level models, random intercepts/slopes, and batch-specific residual variances. Single-study tiers are clearly labelled. Manufacturer tiers contain multiple array designs; sensitivity to finer platform definitions remains untested. |
| Section 7, Plot 1: PCA/UMAP before and after correction | Implemented | Both produced. These use the five-gene panel, not transcriptome-wide QC. |
| Plot 2: five forest plots with study effects and pooled diamonds | Implemented | `pooling_strategy.R:240` onward. |
| Plot 3: ROC/LODO | Excluded from this review | ML excluded by request. |
| Plot 4: expression heatmap across listed indications and healthy controls | Adapted/partial | Current main heatmap shows within-study Hedges' g versus controls, not absolute expression or a separate healthy-control row. AKI is absent from the curated inputs and explicitly disclosed. Separate per-dataset expression heatmaps also exist. |
| Plot 5: expression versus eGFR/fibrosis scatterplots | Missing because targets are unavailable locally | This is a non-ML visualization requirement and remains unmet even with ML excluded. |
| Broad 32-dataset/multi-technology coverage | Reduced scope | Active analysis contains 13 microarray series and selected platforms. For example, the ERCB matrices used contain 71 and 60 samples, rather than the full-series counts of 196 and 195 listed in the annex. The PDF and curated inventory need reconciliation, not an invented 32-study denominator. |

## Findings requiring action

### 1. High: GSE30122 reanalysed controls are treated as independent samples

Locations: `Scripts/stage23_panel_boxplots.R:192`, `:328`, `:411`; `Scripts/pooling_strategy.R:138`; `Scripts/stage4_MF_stat_mod.R:300`; `Scripts/mega_analysis.R:139`; `Scripts/summary_figures.R:177`.

The series matrix explicitly links **25 control reanalyses** to original GSMs: 13 glomerular and 12 tubular. For example, GSM758485 is a reanalysis of GSM757001. Both enter the pipeline as separate samples because the handoff retains GSM and study but no canonical specimen identifier or reanalysis relationship.

GSE30122 therefore has 26 glomerular control records representing 13 linked original/reanalysis pairs, and 24 tubular control records representing 12 pairs. The values differ: the largest panel-gene difference between linked records is approximately 0.493 in glomeruli and 2.298 in tubules. They are alternative analyses of linked source samples, not additional independent controls.

Consequences: sample counts overstate independent evidence, within-study SDs mix processing differences with biological variation, effect-size variances use inappropriate independent counts, and Stage 4/mega models assume independent residuals. The reported 850 distinct GSMs must not be interpreted as 850 independent biological samples. Accounting only for these explicit links reduces distinct specimen entities to at most 825, before other replication or patient pairing is considered.

Recommended action: add canonical specimen and patient identifiers, preserve all series for provenance, and define how alternative processed records contribute to each inferential model. Model dependence or specify a canonical representation/sensitivity analysis. The existing mega sensitivity removing two series does not resolve reanalyses inside GSE30122.

### 2. High: known cross-series overlap is still analysed as independent evidence

Locations: `Scripts/pooling_strategy.R:138`, `:177`, `:223`; `Scripts/stage4_MF_stat_mod.R:309`; `Scripts/mega_analysis.R:139`, `:174`.

The 22 GSE30528 GSMs and 22 GSE30529 GSMs also occur in GSE30122. Their five-gene values match exactly. Keeping all 13 series is an explicit project requirement; that does not make these records statistically independent.

The meta-analysis uses independent-effect `rma` models with no sampling covariance. Stage 4 has no specimen-level correlation, and mega-analysis groups random effects by batch, so repeated specimens in different batches remain unlinked. Leave-one-unit-out can leave the omitted unit's samples in another series. Patient pairing across compartments is also not modelled.

A diagnostic rerun of the compartment meta-analysis omitting only the two subset series produced these BH-FDR changes across the five genes within the tubulointerstitial subgroup:

| Gene | Current FDR | Diagnostic FDR without subset series |
| --- | ---: | ---: |
| COL13A1 | 0.0484 | 0.2512 |
| HABP2 | 0.0169 | 0.0986 |
| CFH | 0.0228 | 0.2059 |
| CFI | 0.000766 | 0.00658 |

These are sensitivity results, not corrected estimates: they still retain GSE30122's internal reanalyses. They demonstrate that the overlap affects reported significance materially.

Recommended action: retain all datasets in the inventory and descriptive outputs, but account for dependent evidence in inference; add grouped omission sensitivity that removes an entire overlapping family. The [metafor dependency documentation](https://wviechtb.github.io/metafor/reference/vcalc.html) explains why partial sample overlap violates standard independent-effect models and describes covariance-based approaches.

### 3. High: GSE66494 reference-RNA replicates are counted as eight independent healthy controls

Locations: `Scripts/stage23_panel_boxplots.R:179`, `:233`; `Scripts/pooling_strategy.R:138`, `:426`; `Scripts/stage4_MF_stat_mod.R:235`; `Scripts/mega_analysis.R:130`.

Local descriptions identify discovery controls as normal-kidney total-RNA replicates 1–5 and validation controls as replicates 1–3. The [original study](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0136994) states that commercially purchased normal-kidney RNA supplied the reference. These records do not establish eight independent healthy donors.

The code notices the reference-RNA category in a sensitivity annotation, but the primary analysis still uses `n_Control = 8`, independent-observation models, and an ordinary healthy-control contrast. It also merges discovery and validation sets without a set term. The same reference samples receive inferred individual-sex labels even though reference material may not represent an individual patient.

Consequences: control biological variance and effective sample size are unverified; the resulting g, sampling variance, and healthy-control interpretation are not justified by this design alone. A random-effects model across studies does not repair the within-study reference design.

Recommended action: verify reference donor/lot and replicate identities, represent discovery/validation structure, and report patient-versus-reference findings separately or supply a justified reference-aware model and sensitivity analysis. Merely labelling the study as reference RNA in a plot is insufficient.

### 4. Medium: GSE66494 effects on log2(signal + 1) are labelled ordinary log2 fold changes

Locations: `Scripts/stage23_panel_boxplots.R:306`; `Scripts/plots_per_dataset.R:241`; `Scripts/pooling_strategy.R:149`; `Scripts/mega_analysis.R:244`.

Adding 1 is substantial for median-normalized Agilent values below 1. A difference of means on this transformed scale describes shifted signals, not a fold change in the original signal. Using the same selected probes, diagnostic contrasts are:

| Gene | Current difference on log2(signal + 1) | Difference on log2(signal) |
| --- | ---: | ---: |
| C4BPB | 0.0409 | 0.5208 |
| COL13A1 | 0.2644 | 0.4023 |
| HABP2 | 0.8262 | 1.7043 |

This diagnostic holds probe selection fixed and is not a proposed replacement preprocessing pipeline. It demonstrates sensitivity to the arbitrary offset. Hedges' g is invariant to linear rescaling, not to this nonlinear transformation.

Recommended action: justify the low-signal transformation and its offset, validate sensitivity, and label the resulting estimand accurately before pooling it with ordinary log2-expression differences.

### 5. Medium: standalone GSE30122 plots test across tissue compartments

Locations: `Scripts/plots_per_dataset.R:151`, `:260`, `:392`.

The registry sets `Ctrl-Glom` as a universal reference; the appropriate within-compartment `pairs` setting is commented out. Consequently the plotted tests include Ctrl-Glom versus Ctrl-Tub and Ctrl-Glom versus DKD-Tub. The latter combines disease and compartment differences and is not a tubular disease-versus-control test.

Recommended action: use the explicit glomerular and tubular comparison pairs. Also identify these standalone Wilcoxon p-values as unadjusted exploratory results or apply an explicit multiple-testing policy. These plots are separate from the main compartment-aware meta-analysis.

### 6. Medium: agreement plots overstate validation as independent

Location: `Scripts/summary_figures.R:35`, `:137`.

Stage 4 and DL pooling are fitted to the same expression observations and share overlap, labels, and preprocessing. Their correlation is useful method concordance, but it is not independent validation. The per-study comparison also contrasts sex-adjusted and unadjusted quantities. The Stage 4 per-study effects exported at `stage4_MF_stat_mod.R:345` are not inputs to the pooling script, despite the preceding comment saying they will be pooled.

Recommended action: describe the figure as concordance between modelling approaches on shared data, explicitly identify the adjustment differences, and correct the misleading handoff comment.

## Correct implementation choices worth preserving

- MaxMean probe collapse uses the requested WGCNA function and all five genes are retained.
- Positional expression/metadata alignment avoids accidental first-match indexing of duplicate GSMs.
- The batch-correction design is checked for full rank and preserves disease effects by compartment. Omitting an aliased compartment main effect here is reasonable, rather than an implementation error. Cross-compartment absolute levels are not preserved.
- Stage 4 models batch explicitly using pre-correction values. This follows the [limma guidance](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf) to include batch factors in inferential models rather than treating corrected values as original observations.
- Primary standardized-effect pooling is useful for heterogeneous study scales; the SMD calculation, bias correction, and DL method match the plan, conditional on valid independent observations and comparison groups.
- Mega-analysis includes random disease slopes and heterogeneous residual variances. Every current multi-study fit converged with the random-slope model; no intercept-only fallback was used.
- Missing clinical outcomes are disclosed rather than fabricated. Sex inference has observed concordance of 40/40 and 60/61 in the two cohorts where both usable inference and reported sex exist; this supports those checks but does not validate every platform or replace missing clinical endpoints.

## Additional limitations and robustness improvements

- `stage23_panel_boxplots.R:429` computes SD across both conditions, not a fitted residual SD. Its filenames and comments call it residual spread, but disease effects and case/control balance contribute to this statistic. Stage 4 independently estimates actual model-residual variances.
- Stage 4 estimates one common disease coefficient per compartment with plug-in variance weights and ordinary WLS covariance. It does not accommodate study-specific disease heterogeneity as the mega-analysis does. This is a modelling limitation, not a direct violation of a specified PDF estimator.
- Missing-file or annotation failures can skip a cohort; missing genes can be dropped by intersection. Only `length(panels) > 1` is asserted. The current run includes all 13 × five genes, but production checks should enforce that invariant explicitly.
- Batch R² reductions support removal of location differences in this five-gene representation. PCA/UMAP mixing alone does not demonstrate preservation of all biological signal or validate downstream p-values.
- A 16-series curated selection is internally consistent with the current 13 included plus three excluded series. The external inventory and full search provenance were not available for verification, so full PRISMA/search compliance is not established.

## Verification performed

1. Extracted the complete PDF text and visually inspected the pages specifying Stages 2–4, pooling, and visualizations; consulted `CLAUDE.md` and `AGENTS.md` for authorized scope.
2. Parsed all eight R scripts successfully under R 4.5.2. Statically inspected the active pipeline and standalone plotting script. The retired Stage 1 and network diagnostic were not executed.
3. Executed the five active pipeline scripts, redirecting only configuration output paths to `/tmp/renorisk_audit_results`; input paths were pointed at this checkout. No statistical code was patched.
4. All five scripts completed successfully. All **27 regenerated workbooks matched existing workbook contents numerically** at tolerance `1e-10`. This establishes reproducibility of these outputs, not statistical correctness.
5. Confirmed 894 observation records per gene, 13 series, 14 batch units, five genes, and no missing panel expression values. Parsed local sample relationships to identify the additional 25 reanalyses and confirmed exact expression equality for the 44 cross-series repeated GSMs.
6. Ran temporary diagnostic meta-analysis and transformation calculations reported above. Existing results were left intact. Expected large-SMD warnings and an S4Vectors R-build-version warning occurred; neither stopped execution.

## Priority order

First resolve biological specimen identity, reanalysis relationships, reference-RNA replication, and inferential dependence while retaining the required series inventory. Then validate transformation choices and fix standalone comparisons. Finally update the protocol compliance statement to explicitly document retired Stage 1, reduced cohort/platform coverage, unavailable clinical covariates and scatterplots, and the adapted indication heatmap. ML remains excluded from this audit.
