# Mitochondrial respiratory capacity and coupling dynamics in sepsis

Analysis code for:

Kopczynski et al. Mitochondrial respiratory capacity and coupling dynamics discriminate mortality risk in patients with sepsis. *Science Advances* (submitted).

The code reproduces every table, figure and statistic reported in the article and its Supplementary Materials from the study data set (40 patients with sepsis requiring vasopressor support; 37 with paired respirometry at admission and 72 h).

## Data

The data file is not distributed with the code; access is described in the Data and materials availability statement of the article. Place it in `data/` as `GDF15final.csv` (comma- or semicolon-delimited) or `GDF15final.xlsx`. Variables are listed in `data/README.md`.

## Requirements

R 4.3 or later and the packages

```r
install.packages(c("dplyr", "tidyr", "tibble", "ggplot2", "patchwork", "ggrepel", "scales",
                   "pROC", "logistf", "cluster", "mclust", "fpc", "diptest",
                   "readxl", "stringr", "officer", "flextable"))
```

## Running

From the repository root:

```
Rscript run_all.R
```

or run the scripts in `R/` individually, in order, after `R/01_prepare_data.R`. Results are written to `output/`, one folder per script. The full run takes several minutes, most of it in the repeated cross-validation.

## Scripts and outputs

| Script | Article | Content |
|---|---|---|
| `00_setup.R` | | Packages, paths, data loading, derived variables, shared figure style |
| `01_prepare_data.R` | | Reads the raw file, standardises names and types, writes `data/analysis_data.rds` |
| `02_table1.R` | Table 1 | Baseline characteristics; counts quoted in Results |
| `03_figures_1_3.R` | Figs. 1 to 3 | Mediators and respirometry by 28-day and 6-month mortality, Mann-Whitney U tests; admission IL-10 and IL-6 medians |
| `04_figure4_figS1.R` | Fig. 4, fig. S1 | Principal components and marker projection (4A); cross-validated ROC and AUC at 28 days and 6 months (4B to 4E); mortality by baseline coupling and direction of change with Fisher exact tests (4F); coupling trajectories (fig. S1) |
| `05_dimensions_recovery.R` | Results | Correlation of Complex I with coupling and with Complex II; recovery of coupling predicted by admission Complex I (AUC, Spearman, Firth model with SOFA); Complex I coefficient adjusted for the coupling trajectory |
| `06_mitokine_correlations.R` | Results | 32 Spearman correlations of GDF15, FGF21, leptin and ghrelin with respirometry at 0 h and 72 h, FDR corrected |
| `07_subgroup_tests.R` | Results | Gap statistic, silhouette width, Gaussian mixture model, bootstrap stability and dip tests |
| `08_table2_supplementary_tables.R` | Table 2, tables S1 to S3 | Cross-validated AUC, comparisons with SOFA, regression to the mean, ANCOVA, recovery of coupling, 37 versus 40 patients |
| `09_table3_nested_tests.R` | Table 3 | Penalised likelihood ratio tests and cross-validated gain in AUC |
| `10_figure5_counterfactual.R` | Fig. 5 | Counterfactual trial-enrichment analysis |

## Methods implemented

- Mortality horizons: ICU (`mort_uti`), in-hospital (`mort_hosp`), 28 days (`mort_28d`) and 6 months (`mort_6m`).
- Coupling efficiency (BCE) = (P − Ln)/P. Coupling trajectory = 100 × (BCE 72 h − BCE 0 h) / BCE 0 h.
- Discrimination: AUC from five-fold cross-validation repeated 20 times, with Firth penalised logistic regression (`logistf`, `pl = FALSE`) fitted inside each training fold; out-of-fold predictions are pooled within each repetition and the AUC is averaged across repetitions. Apparent AUCs are not computed because several models separate completely.
- Nested comparisons: `logistftest` with the added coefficients constrained to zero within the full Firth model, standardised predictors.
- Paired AUC comparisons: DeLong test on out-of-fold predicted probabilities.
- Counterfactual analysis: single-variable logistic models, absolute risk reduction and number needed to treat, 3,000 bootstrap resamples with refitting.

## Reproducibility

Each script sets `set.seed(42)`. Cross-validation folds and bootstrap resamples depend on the order in which models and horizons are processed; that order is fixed in the scripts and should not be changed if the published values are to be reproduced exactly. `output/sessionInfo.txt` records the package versions used in a run.

## License

MIT (see `LICENSE`).
