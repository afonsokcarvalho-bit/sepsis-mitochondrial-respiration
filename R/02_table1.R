# Table 1. Baseline characteristics of the study population (n = 40).
# Mean (SD) or median (interquartile range) as in the published table;
# categorical variables as n (%). Quartiles use the R default (type 7).

source("R/00_setup.R")
df <- load_data()
OUT <- out_dir("table1")

# Columns for renal replacement therapy and source of infection.
COL_RRT    <- "hd"            # renal replacement therapy, 0/1
COL_SOURCE <- NA_character_   # source of infection (not in the data file)

ROWS <- tibble::tribble(
  ~var,             ~label,                                              ~type,
  "idade",          "Age, years",                                        "mean",
  "sex_masc",       "Sex (male)",                                        "n",
  "bmi",            "Body mass index, kg/m2",                            "mean",
  "saps_3",         "SAPS 3 score",                                      "mean",
  "sofa_0h",        "SOFA score, 0 h",                                   "median",
  "sofa_72h",       "SOFA score, 72 h",                                  "median",
  "fc_max_0h",      "Maximum heart rate, 0 h, beats/min",                "mean",
  "ench_cap_0h",    "Capillary refill time, 0 h, s",                     "median",
  "bh_cum_0h",      "Cumulative fluid balance, 0 h, mL",                 "median",
  "vasopr_max_0h",  "Maximum norepinephrine dose, 0 h, ug/kg/min",       "median",
  "lactate_0h",     "Arterial lactate at diagnosis, mmol/L",             "median",
  "svo2_0h",        "Central venous oxygen saturation, %",               "median",
  "gapco2_0h",      "p(v-a)CO2 at diagnosis, mmHg",                      "median",
  "troponine_0h",   "Troponin, 0 h, ng/L",                               "median",
  "albumin_0h",     "Albumin, 0 h, g/dL",                                "median",
  "pf_0h",          "PaO2/FiO2 ratio, 0 h",                              "mean",
  "crp_0h",         "C-reactive protein, 0 h, mg/L",                     "median",
  "crp_72h",        "C-reactive protein, 72 h, mg/L",                    "median",
  "leuc_0h",        "Total leukocyte count, 0 h, cells/mm3",             "median",
  "leuc_72h",       "Total leukocyte count, 72 h, cells/mm3",            "median",
  "vm",             "Mechanical ventilation",                            "n",
  COL_RRT,          "Renal replacement therapy",                         "n",
  "corticoide_0h",  "Corticosteroid therapy, 0 h",                       "n",
  "mort_uti",       "Mortality, ICU",                                    "n",
  "mort_hosp",      "Mortality, in-hospital",                            "n",
  "mort_28d",       "Mortality, 28 days",                                "n",
  "mort_6m",        "Mortality, 6 months",                               "n")

summ <- function(x, type) {
  x <- suppressWarnings(as.numeric(as.character(x))); n <- sum(!is.na(x))
  dg <- if (stats::median(x, na.rm = TRUE) >= 100) 0 else 1
  val <- switch(type,
    mean   = sprintf("%s (%s)", fmt(mean(x, na.rm = TRUE), dg), fmt(sd(x, na.rm = TRUE), dg)),
    median = { q <- stats::quantile(x, c(.5, .25, .75), na.rm = TRUE)
               sprintf("%s (%s to %s)", fmt(q[1], dg), fmt(q[2], dg), fmt(q[3], dg)) },
    n      = sprintf("%d (%.1f%%)", sum(x == 1, na.rm = TRUE), 100 * mean(x == 1, na.rm = TRUE)))
  c(val, n)
}

T1 <- dplyr::bind_rows(lapply(seq_len(nrow(ROWS)), function(i) {
  v <- ROWS$var[i]
  if (is.na(v) || !v %in% names(df))
    return(tibble::tibble(Characteristic = ROWS$label[i], Value = "column not set", n = NA_integer_))
  s <- summ(df[[v]], ROWS$type[i])
  tibble::tibble(Characteristic = ROWS$label[i], Value = s[1], n = as.integer(s[2]))
}))

if (!is.na(COL_SOURCE) && COL_SOURCE %in% names(df)) {
  src <- as.data.frame(table(df[[COL_SOURCE]], useNA = "ifany"))
  T1 <- dplyr::bind_rows(T1, tibble::tibble(
    Characteristic = paste0("Source of infection: ", src$Var1),
    Value = sprintf("%d (%.1f%%)", src$Freq, 100 * src$Freq / nrow(df)),
    n = nrow(df)))
}
write_out(T1, OUT, "table1.csv")

# Values quoted in the Results
lac <- df$lactate_0h[!is.na(df$lactate_0h)]
txt <- c(
  sprintf("Patients: %d", nrow(df)),
  sprintf("Arterial lactate < 2 mmol/L at inclusion: %d of %d (%.0f%%)",
          sum(lac < 2), length(lac), 100 * mean(lac < 2)),
  sprintf("Paired respirometry (0 h and 72 h): %d", sum(paired(df))),
  "",
  "Observations available, where fewer than all patients:",
  sprintf("  %-16s %d", T1$Characteristic[!is.na(T1$n) & T1$n < nrow(df)],
          T1$n[!is.na(T1$n) & T1$n < nrow(df)]),
  sprintf("  %-16s %d", c("lactate_24h", "lactate_72h"),
          c(sum(!is.na(df$lactate_24h)), sum(!is.na(df$lactate_72h)))))
writeLines(txt, file.path(OUT, "table1_notes.txt"))
print(as.data.frame(T1), row.names = FALSE)
