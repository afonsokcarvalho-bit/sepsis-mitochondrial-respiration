# Figures 1 to 3. Inflammatory mediators (Fig. 1), perfusion, endocrine and
# metabolic markers (Fig. 2) and mitochondrial respiration (Fig. 3) at
# admission and 72 h (lactate: 24 h), by 28-day and 6-month mortality.
# Survivors and non-survivors compared with the Mann-Whitney U test (exact
# P values when there are no ties, normal approximation otherwise).
# Output: output/figures_1_3/ (figures and figures_1_3_statistics.csv)

source("R/00_setup.R")
df <- load_data()
OUT <- out_dir("figures_1_3")
set.seed(42)   # jitter only

df <- df %>% mutate(
  Outcome_28d   = factor(ifelse(mort_28d == 0, "Surv.", "Non-Surv."), levels = c("Surv.", "Non-Surv.")),
  Outcome_6m    = factor(ifelse(mort_6m  == 0, "Surv.", "Non-Surv."), levels = c("Surv.", "Non-Surv.")),
  tnfa_il10_0h  = tnfa_0h / il10_0h,     tnfa_il10_72h = tnfa_72h / il10_72h,
  lep_ghr_0h    = leptin_0h / ghrelin_0h, lep_ghr_72h  = leptin_72h / ghrelin_72h)
for (s in c("_0h", "_72h")) {               # IL-1 beta is stored as il1_*
  if (paste0("il1", s) %in% names(df)) df[[paste0("il1b", s)]] <- df[[paste0("il1", s)]]
}

stars <- function(p) dplyr::case_when(p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", TRUE ~ "ns")

panel <- function(var_t1, var_t2, outcome, y_label, title, time2 = "72h") {
  is_28d <- outcome == "Outcome_28d"
  fill <- if (is_28d) "#2CA02C" else "#FF7F0E"
  d <- df %>%
    dplyr::select(Outcome = !!sym(outcome), T1 = !!sym(var_t1), T2 = !!sym(var_t2)) %>%
    pivot_longer(c(T1, T2), names_to = "Time", values_to = "Value") %>%
    mutate(Time = factor(ifelse(Time == "T1", "0h", time2), levels = c("0h", time2))) %>%
    drop_na(Value, Outcome)
  pv <- d %>% group_by(Time) %>%
    summarise(p = suppressWarnings(wilcox.test(Value ~ Outcome)$p.value), .groups = "drop") %>%
    mutate(sig = stars(p)) %>% filter(sig != "ns")
  p <- ggplot(d, aes(Outcome, Value, fill = Outcome)) +
    geom_violin(alpha = 0.3, colour = NA, trim = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.8, colour = "black", outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, colour = "black", shape = 21) +
    scale_fill_manual(values = c("Surv." = fill, "Non-Surv." = fill)) +
    facet_wrap(~ Time, strip.position = "bottom") +
    labs(title = title, y = y_label, x = NULL) +
    theme_classic(base_family = "serif", base_size = 14) +
    theme(plot.tag = element_text(face = "bold", size = 20),
          plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
          strip.background = element_blank(), strip.placement = "outside",
          strip.text = element_text(face = "bold", size = 14, colour = "black"),
          axis.text = element_text(colour = "black", size = 12),
          axis.title.y = element_text(face = "bold", size = 13),
          legend.position = "none")
  if (nrow(pv)) {
    top <- max(d$Value); rng <- top - min(d$Value)
    yb <- top + rng * 0.05; tick <- rng * 0.02
    p <- p +
      geom_segment(data = pv, aes(x = 1, xend = 2, y = yb, yend = yb), inherit.aes = FALSE) +
      geom_segment(data = pv, aes(x = 1, xend = 1, y = yb, yend = yb - tick), inherit.aes = FALSE) +
      geom_segment(data = pv, aes(x = 2, xend = 2, y = yb, yend = yb - tick), inherit.aes = FALSE) +
      geom_text(data = pv, aes(x = 1.5, y = top + rng * 0.07, label = sig),
                inherit.aes = FALSE, size = 6, fontface = "bold") +
      coord_cartesian(ylim = c(min(d$Value), top + rng * 0.15))
  }
  p
}

FLUX <- "O\u2082 flux (pmol\u00b7s\u207b\u00b9\u00b7mg\u207b\u00b9)"
SPEC <- tibble::tribble(
  ~fig, ~v0,            ~v72,            ~ylab,                   ~name,           ~t2,
  1,    "tnfa_0h",      "tnfa_72h",      "TNF-\u03b1 (pg/mL)",    "TNF-\u03b1",    "72h",
  1,    "il1b_0h",      "il1b_72h",      "IL-1\u03b2 (pg/mL)",    "IL-1\u03b2",    "72h",
  1,    "il6_0h",       "il6_72h",       "IL-6 (pg/mL)",          "IL-6",          "72h",
  1,    "crp_0h",       "crp_72h",       "CRP (mg/L)",            "CRP",           "72h",
  1,    "il10_0h",      "il10_72h",      "IL-10 (pg/mL)",         "IL-10",         "72h",
  1,    "tnfa_il10_0h", "tnfa_il10_72h", "TNF-\u03b1/IL-10 ratio", "TNF-\u03b1/IL-10", "72h",
  2,    "lactate_0h",   "lactate_24h",   "Lactate (mmol/L)",      "Lactate",       "24h",
  2,    "leptin_0h",    "leptin_72h",    "Leptin (pg/mL)",        "Leptin",        "72h",
  2,    "ghrelin_0h",   "ghrelin_72h",   "Ghrelin (pg/mL)",       "Ghrelin",       "72h",
  2,    "gdf15_0h",     "gdf15_72h",     "GDF15 (pg/mL)",         "GDF15",         "72h",
  2,    "fgf21_0h",     "fgf21_72h",     "FGF21 (pg/mL)",         "FGF21",         "72h",
  2,    "lep_ghr_0h",   "lep_ghr_72h",   "Leptin/ghrelin ratio",  "Leptin/ghrelin", "72h",
  3,    "r_0h",         "r_72h",         FLUX,                    "LEAK",          "72h",
  3,    "ci_0h",        "ci_72h",        FLUX,                    "Complex I",     "72h",
  3,    "cii_0h",       "cii_72h",       FLUX,                    "Complex II",    "72h",
  3,    "bce_0h",       "bce_72h",       "(P \u2212 Ln)/P",       "Coupling",      "72h")

# each marker gives two panels: 28 days (left) and 6 months (right)
PANELS <- list(); STATS <- list()
for (i in seq_len(nrow(SPEC))) {
  s <- SPEC[i, ]
  for (oc in c("Outcome_28d", "Outcome_6m")) {
    hz <- if (oc == "Outcome_28d") "28 DAYS" else "6 MONTHS"
    PANELS[[length(PANELS) + 1]] <- list(fig = s$fig,
      p = panel(s$v0, s$v72, oc, s$ylab, paste(toupper(s$name), hz), s$t2))
    for (tp in c(s$v0, s$v72)) {
      x <- df[[tp]]; g <- df[[oc]]; ok <- !is.na(x) & !is.na(g)
      q <- function(k) { z <- x[ok & g == k]
        sprintf("%.2f (%.2f to %.2f)", median(z), quantile(z, .25), quantile(z, .75)) }
      STATS[[length(STATS) + 1]] <- tibble::tibble(
        figure = s$fig, marker = s$name, horizon = hz,
        time = if (tp == s$v0) "0 h" else sub("h", " h", s$t2),
        n = sum(ok), survivors = q("Surv."), non_survivors = q("Non-Surv."),
        P = suppressWarnings(wilcox.test(x[ok] ~ g[ok])$p.value))
    }
  }
}
STATS <- dplyr::bind_rows(STATS) %>% mutate(sig = stars(P))
write_out(STATS, OUT, "figures_1_3_statistics.csv")

assemble <- function(k, ncol = 4) {
  ps <- lapply(Filter(function(z) z$fig == k, PANELS), `[[`, "p")
  patchwork::wrap_plots(ps, ncol = ncol) + patchwork::plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 20))
}
ggsave(file.path(OUT, "Figure1.png"), assemble(1), width = 16, height = 12, dpi = 300)
ggsave(file.path(OUT, "Figure2.png"), assemble(2), width = 16, height = 12, dpi = 300)
ggsave(file.path(OUT, "Figure3.png"), assemble(3), width = 16, height = 8,  dpi = 300)

# Admission IL-10 and IL-6, all patients (Results, inflammatory response)
writeLines(sprintf("Admission medians, all patients: IL-10 %.1f pg/mL, IL-6 %.1f pg/mL",
                   median(df$il10_0h, na.rm = TRUE), median(df$il6_0h, na.rm = TRUE)),
           file.path(OUT, "il10_il6_admission.txt"))
