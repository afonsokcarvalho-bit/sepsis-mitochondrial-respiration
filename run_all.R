# Runs the complete analysis from the repository root, each script in a
# fresh R session:  Rscript run_all.R

scripts <- c("R/01_prepare_data.R",
             "R/02_table1.R",
             "R/03_figures_1_3.R",
             "R/04_figure4_figS1.R",
             "R/05_dimensions_recovery.R",
             "R/06_mitokine_correlations.R",
             "R/07_subgroup_tests.R",
             "R/08_table2_supplementary_tables.R",
             "R/09_table3_nested_tests.R",
             "R/10_figure5_counterfactual.R")

rscript <- file.path(R.home("bin"), "Rscript")
for (s in scripts) {
  message("== ", s)
  if (system2(rscript, s) != 0) stop("Failed: ", s)
}

dir.create("output", showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")
