# Full-expression sample-overlap audit

Date: 2026-10-19; re-run 2026-11-04 over the 19 series now in use. This is an
exploratory provenance audit. No samples were removed, labels changed, classifier
trained, or existing analysis outputs overwritten.

## Re-run 2026-11-04

All 19 series in use (14 microarray + 5 RNA-Seq, 1322 samples; the `-dup` and `-excl`
folders are not read). The 8 series added since 2026-10-19 (GSE93798, GSE35488,
GSE111154, GSE142025, GSE141295, GSE175759, GSE137570, GSE219185) produce **no
candidate pair**: their best reciprocal-nearest cross-study correlation is 0.73
(GSE142025 / GSE137570), far below the 0.95 screen. The findings below are unchanged:
the 100 candidates are the same 96 GSE108109 / GSE200828 pairs and 4 GSE66494
reference-RNA pairs; GSE121233 / GSE129973 still show 10 reciprocal-nearest pairs at
0.854–0.898 (below the screen). Known overlaps: 25 documented GEO relations within
GSE30122, none above the screen. GSE30528/GSE30529, the former positive controls, are
no longer read (dropped 2026-11-03). Partial overlaps are kept as they are (user
instruction, 2026-11-03).

## Main findings

- **GSE108109 versus GSE200828: 96 strong candidate reused-sample pairs.** All
  pairs are reciprocal nearest neighbours across 25,582 mapped genes. Correlations
  of per-study, per-gene standardized profiles range from 0.9682101 to 0.9901481;
  raw gene-profile Spearman correlations range from 0.9967175 to 0.9980959.
  Both disjoint gene panels exceed the screening threshold for all 96 pairs.
  All pairs share the numeric identifier embedded in their sample titles;
  diagnosis text agrees for 94/96. These are different GSM accessions, and their
  deposited native-probe values are not exactly identical. The evidence supports
  sample reuse/reprocessing, not a claim of byte-identical files or independently
  verified patient identity. The local series summaries identify NEPTUNE in both.
- **GSE121233 versus GSE129973: 10 additional candidate overlaps requiring review.**
  Every GSE121233 sample has a different reciprocal-nearest counterpart in
  GSE129973. Standardized correlations range from 0.8542724 to 0.8975308, with
  consistent correlations in both disjoint gene panels (0.8563506–0.8969782 and
  0.8511353–0.8980410). These do not exceed the primary 0.95 retrieval threshold;
  find them in `nearest_and_screened`, filtering these two studies and
  `reciprocal_nearest = TRUE`. This is suggestive, not confirmed identity.
- Four additional pairs exceeding the primary screen occur among GSE66494
  control/reference-RNA profiles. Do not interpret these as four independent
  patient duplicates: the pairs involve five GSM records, and the normal-kidney
  total-RNA source/replicate structure needs investigation.

Examples of the first finding:

| GSE108109 | GSE200828 | Shared title identifier |
| --- | --- | --- |
| GSM2889867: Living donor [B1_Glom_LD_3] | GSM6044637: Glomerular_3_Living donor | 3, living donor |
| GSM2889863: FSGS [B1_Glom_FSGS_87] | GSM6044614: Glomerular_87_FSGS | 87, FSGS |
| GSM2889872: MCD [B1_Glom_MCD_128] | GSM6044595: Glomerular_128_MCD | 128, MCD |

Diagnosis discrepancies among these 96 strong candidate pairs:

| GSE108109 | GSE200828 |
| --- | --- |
| GSM2889851: FSGS, identifier 4 | GSM6044606: MCD, identifier 4 |
| GSM2889875: MCD, identifier 95 | GSM6044611: FSGS, identifier 95 |

The numeric suffix alone is not a globally unique patient identifier. Diagnosis,
cohort provenance and whole-expression matching provide additional evidence.

## Methods and reproducibility

Run from the project root:

```sh
Rscript Scripts/by-claude/audit_expression_overlap.R
```

The script reads only `Data/Raw_filtered/`: the 19 series in use, 1322 samples
(RNA-Seq values from `<id>_expression.txt`; symbols / Ensembl IDs mapped to Entrez
with org.Hs.eg.db; GSE35488's bare Entrez IDs accepted). It evaluates 873,181
sample pairs in 190 within/between-study comparisons, using all shared mapped
genes (9,007–32,920, depending on the study pair). (2026-10-19 run: 13 matrices,
894 entries incl. GSE30528/GSE30529 as positive controls, 399,171 pairs.)

Map native features to Entrez IDs using installed annotation databases. Exclude
ambiguous mappings and nonfinite rows. Collapse multiple probes by their mean
for the fingerprint audit only; the RenoRisk panel's MaxMean extraction remains
unchanged. Transform GSE66494 with the existing log2(x + 1) convention.

Calculate Pearson correlations after gene centering within each study and after
gene centering plus gene-SD scaling within each study. Calculate raw-profile
Spearman correlations. Repeat standardized correlations on two disjoint gene
sets. Retrieve pairs with standardized r >= 0.95 in the full set and both gene
sets. This is a broad heuristic screen, not a calibrated identity decision rule.
Also retain nearest neighbours, known identical GSMs and explicit GEO sample
relations. For selected pairs with matching probe namespaces, compare deposited
native-probe values directly (equal fraction, RMSE and maximum absolute error).
No disease labels enter the expression calculations. This audit standardization
does not define a classifier preprocessing procedure.

Workbook sheets:

- `inventory`: input dimensions and mapping coverage.
- `pair_summary`: all study-pair comparisons and correlation background summaries.
- `new_candidates`: 100 undocumented high-similarity pairs among the 1322 GSMs:
  96 NEPTUNE-related pairs and four reference-RNA pairs.
- `known_overlaps`: identical GSMs and explicit GEO relations (25, within GSE30122).
- `nearest_and_screened`: all retrieved pairs, including lower-scoring neighbours.

A targeted check of the 10 HTA candidate pairs additionally re-centred the
70,523 shared native probes within the matched ten-sample subsets. Matched
correlations were 0.9176154–0.9322143, versus a maximum off-diagonal correlation
of 0.1388853. None of their deposited probe values were exactly equal. Since the
matched subset was selected using expression similarity, these additional
numbers are descriptive follow-up, not an independent validation or p-value.

## Sensitivity and interpretation

The screen recovers all 44 known identical-GSM cross-series matches at r >= 0.95.
However, the 25 explicitly linked reanalysis pairs inside GSE30122 do not pass
that stringent threshold (r ranges from 0.03210546 to 0.9454461). Therefore the
screen cannot establish that lower-scoring pairs are independent. Gene-panel
agreement is a stability check, not an independent patient-identity assay.

Submitter-processed microarrays differ in platform, normalization, gene centering,
and tissue compartment. Whole-expression comparisons can detect profile reuse,
but can miss the same donor sampled in different compartments or times. Similar
control reference RNA can also resemble technical replication without identifying
individual donors. Confirmation needs specimen/patient crosswalks, original-array
identifiers when available, or clarification from the original investigators.

## Consequence for classifier validation

850 remains the count of unique GSM accessions, not independent specimens or
patients. Do not publish a revised independent-patient count by subtracting
candidate-pair counts. Retain all data during review. For conservative LOGO
validation, keep GSE108109 and GSE200828 in one group; also evaluate grouping
GSE121233 with GSE129973 while the weaker overlap is investigated. Fold checks
must use resolved specimen/patient groups, not GSM uniqueness alone. These are
recommendations; this audit does not modify the pipeline or its grouping rules.

Methodological background: Waldron et al., *The Doppelgänger Effect: Hidden
Duplicates in Databases of Transcriptome Profiles*, JNCI (2016),
https://doi.org/10.1093/jnci/djw146. This script is a transparent exploratory
correlation screen, not an implementation of the doppelgangR algorithm.
