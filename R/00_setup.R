# Shared settings for all analysis scripts.
# Scripts are run from the repository root, e.g. source("R/04_figure4_figS1.R").

PKGS <- c("dplyr", "tidyr", "tibble", "ggplot2", "patchwork", "ggrepel", "scales",
          "pROC", "logistf", "cluster", "mclust", "fpc", "diptest",
          "readxl", "stringr", "officer", "flextable")
missing_pkgs <- PKGS[!vapply(PKGS, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs))
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall with install.packages(c(\"",
       paste(missing_pkgs, collapse = "\", \""), "\"))")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})
options(OutDec = ".", scipen = 6)

DIR_DATA <- "data"
DIR_OUT  <- "output"
F_RAW    <- file.path(DIR_DATA, "GDF15final.csv")
F_RAW_XL <- file.path(DIR_DATA, "GDF15final.xlsx")
F_DATA   <- file.path(DIR_DATA, "analysis_data.rds")

out_dir <- function(name) {
  d <- file.path(DIR_OUT, name)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

# Mortality outcomes. Several scripts process them in increasing order of
# events (ICU, 28 days, in-hospital, 6 months); that order fixes the sequence of
# random numbers and therefore the published cross-validated values.
OUTCOMES <- c(ICU = "mort_uti", `In-hospital` = "mort_hosp",
              `28 days` = "mort_28d", `6 months` = "mort_6m")
HORIZONS <- names(OUTCOMES)                     # display order
nice <- function(v) names(OUTCOMES)[match(v, OUTCOMES)]
by_events <- function(df) OUTCOMES[order(vapply(OUTCOMES,
  function(v) sum(df[[v]], na.rm = TRUE), numeric(1)))]

# Analysis data set with derived variables.
#   d_<p>    absolute change from admission (0 h) to 72 h
#   rel_<p>  relative change, 100 * (72 h - 0 h) / 0 h; rel_bce is the
#            coupling trajectory
#   d_il10   change in IL-10 from 0 h to 72 h
#   lepghre_0h  leptin to ghrelin ratio at admission
load_data <- function() {
  if (!file.exists(F_DATA)) stop("Run R/01_prepare_data.R first.")
  df <- readRDS(F_DATA)
  names(df) <- trimws(names(df))
  if ("age" %in% names(df) && !"idade" %in% names(df)) df <- dplyr::rename(df, idade = age)
  if ("lactate_24h_0h" %in% names(df) && !"lactate_24h" %in% names(df))
    df <- dplyr::rename(df, lactate_24h = lactate_24h_0h)
  for (v in OUTCOMES) df[[v]] <- as.numeric(as.character(df[[v]]))
  for (p in c("r", "ci", "cii", "bce")) {
    a <- df[[paste0(p, "_0h")]]; b <- df[[paste0(p, "_72h")]]
    df[[paste0("d_", p)]]   <- b - a
    df[[paste0("rel_", p)]] <- ifelse(a > 0, 100 * (b - a) / a, NA)
  }
  df$d_il10     <- df$il10_72h - df$il10_0h
  df$lepghre_0h <- ifelse(df$ghrelin_0h > 0, df$leptin_0h / df$ghrelin_0h, NA)
  df$.row <- seq_len(nrow(df))
  df
}

# Patients with paired respirometry (alive and sampled at 72 h)
paired <- function(df) !is.na(df$bce_0h) & !is.na(df$bce_72h)

# Wilson 95% interval, in percent
wilson <- function(x, n) {
  if (n == 0) return(c(NA, NA))
  p <- x / n; z <- 1.96; d <- 1 + z^2 / n
  c(100 * max((p + z^2/(2*n) - z * sqrt(p*(1-p)/n + z^2/(4*n^2))) / d, 0),
    100 * min((p + z^2/(2*n) + z * sqrt(p*(1-p)/n + z^2/(4*n^2))) / d, 1))
}

fmt <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = d))
fmt_p <- function(x) ifelse(is.na(x), "NA",
  ifelse(x < 1e-4, formatC(x, format = "e", digits = 1),
  ifelse(x < 0.001, formatC(x, format = "f", digits = 4),
         formatC(x, format = "f", digits = 3))))

write_out <- function(d, dir, file)
  write.csv(d, file.path(dir, file), row.names = FALSE, fileEncoding = "UTF-8")

# Shared figure style (Figures 4 and 5)
FONT  <- "serif"
GREEN <- "#49AD49"; ORANG <- "#FF9130"; TEAL <- "#2E8B8B"; PLUM <- "#8E6C9B"
SAND  <- "#C9A227"; SAND_L <- "#EFDFA8"; BROWN <- "#B5651D"
BLACK <- "#1A1A1A"; GREY <- "#9A9A9A"; REF <- "#2B2B2B"

theme_paper <- function()
  ggplot2::theme_classic(base_size = 12, base_family = FONT) +
  ggplot2::theme(
    axis.line  = ggplot2::element_line(colour = BLACK, linewidth = .45),
    axis.ticks = ggplot2::element_line(colour = BLACK, linewidth = .45),
    axis.text  = ggplot2::element_text(colour = BLACK, size = 10, family = FONT),
    axis.title = ggplot2::element_text(colour = BLACK, size = 11, family = FONT),
    plot.title = ggplot2::element_text(size = 11.5, face = "bold", hjust = .5,
                                       family = FONT, colour = BLACK,
                                       margin = ggplot2::margin(b = 6)),
    legend.position = "none",
    plot.margin = ggplot2::margin(7, 9, 5, 7))
TAG <- ggplot2::theme(plot.tag = ggplot2::element_text(size = 20, face = "bold",
         family = FONT, colour = BLACK, hjust = 0, vjust = 1))
LEG <- ggplot2::theme(legend.position = "top",
         legend.title = ggplot2::element_blank(),
         legend.text = ggplot2::element_text(size = 9.5, family = FONT),
         legend.key.size = ggplot2::unit(.38, "cm"),
         legend.margin = ggplot2::margin(0, 0, -4, 0),
         legend.background = ggplot2::element_blank())

save_fig <- function(p, dir, name, w, h) {
  ggplot2::ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h,
                  device = cairo_pdf, bg = "white")
  ggplot2::ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h,
                  dpi = 600, bg = "white")
}
