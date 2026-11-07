# Repository Guidelines

## Project Structure & Module Organization

This R/Bioconductor project validates the RenoRisk five-gene kidney biomarker panel: C4BPA, COL13A1, HABP2, CFH, and CFI.

- `Scripts/`: standalone analysis scripts, organized by pipeline stage.
- `Data/Raw_filtered/<GSE>/`: curated GEO series matrices; the active analysis inputs.
- `Results/`: generated workbooks, figures, and session information, grouped into `Stage23/`, `Stage4/`, `PoolingStrategy/`, `MegaAnalysis/`, and `Summary/`.
- `Examples/`: reference plotting scripts.
- `CLAUDE.md`: detailed analysis decisions and conventions. Consult it alongside the protocol PDF before changing statistical behavior.

## Build, Test, and Development Commands

Run from the project root, in this order:

```bash
Rscript Scripts/stage23_panel_boxplots.R # Map genes and correct batch effects
Rscript Scripts/stage4_MF_stat_mod.R    # Fit multifactor models
Rscript Scripts/pooling_strategy.R     # Pool within-study effects
Rscript Scripts/mega_analysis.R        # Fit platform-tier mixed models
Rscript Scripts/by-claude/summary_figures.R      # Generate agreement and PRISMA summaries
```

Scripts accept no arguments. Configure `cfg$project_dir` for your checkout. Use R 4.5.2 with the required CRAN/Bioconductor packages installed; scripts do not install dependencies. There is no build system.

## Coding Style & Naming Conventions

Follow existing R style: two-space indentation, `<-` assignment, descriptive `snake_case` helpers, and a top-level `cfg` list with derived paths. Preserve numbered sections and comments explaining statistical choices. Use `set.seed(1)` for reproducibility. No formatter or linter is configured.

Write tables as `.xlsx` with `writexl`; read with `readxl`, explicit `col_types`, and `as.data.table()`. Preserve true `NA` values. Keep outputs and `*_sessionInfo.txt` together in each script’s results folder.

## Testing Guidelines

There is no automated test framework or coverage target. Re-run changed stages and downstream consumers. Check `cohort_qc.xlsx` for 13 cohorts × five genes and confirm 894 samples across 14 batch units. Inspect `panel_boxplot_per_cohort.png` for batch alignment. `load_data_test.R` is a GEO connectivity diagnostic, not a unit suite.

## Commit & Pull Request Guidelines

This directory has no Git repository, so no historical commit convention exists. If versioned, use concise imperative subjects. Describe analysis changes, affected outputs, validation commands, and protocol deviations; attach changed figures when relevant. Never commit `Data/Raw/`.

## Analysis Constraints

Use only `Data/Raw_filtered/`; Stage 1 is retired. Preserve all 13 cohorts, including accepted overlaps. Keep one batch-correction model, and use uncorrected expression for Stage 4 and pooling. Document submitter-processing limitations in figures.
