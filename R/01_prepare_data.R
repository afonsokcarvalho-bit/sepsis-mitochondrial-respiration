# Reads the raw data file, standardises column names and types, and saves the
# analysis data set (data/analysis_data.rds) used by every other script.

source("R/00_setup.R")

raw <- tryCatch({
  d <- read.csv(F_RAW, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(d) <= 1)
    d <- read.csv(F_RAW, sep = ";", dec = ",", check.names = FALSE,
                  stringsAsFactors = FALSE)
  d
}, error = function(e) readxl::read_excel(F_RAW_XL))

# column names: spaces to "_", "/" to "_div_", lower case
names(raw) <- tolower(stringr::str_replace_all(names(raw), c(" " = "_", "/" = "_div_")))

# numeric columns stored as text with decimal commas
for (v in setdiff(names(raw), "id"))
  if (is.character(raw[[v]])) raw[[v]] <- as.numeric(stringr::str_replace(raw[[v]], ",", "."))

# patients with known 28-day outcome; outcome coded 0 = survivor, 1 = non-survivor
df <- dplyr::filter(raw, !is.na(mort_28d))
df$mort_28d <- as.numeric(as.factor(df$mort_28d)) - 1
if (max(df$mort_28d, na.rm = TRUE) > 1)
  df$mort_28d <- ifelse(df$mort_28d == min(df$mort_28d, na.rm = TRUE), 0, 1)

saveRDS(df, F_DATA)
message("Saved ", F_DATA, ": ", nrow(df), " patients, ", ncol(df), " variables")
