# Figure 5. Counterfactual trial-enrichment analysis (n = 37, 28-day mortality).
# Single-variable logistic models (Firth only if maximum likelihood separates);
# risk recalculated after a simulated change in the predictor, assuming the
# association is causal.
#   A to C  gains of 10, 20 and 30 percentage points in the coupling trajectory
#   D, E    every target shifted by one SD in the direction that lowers risk
#   F       30-point gain applied to all patients or only to those whose
#           coupling did not increase
# Intervals: 3,000 bootstrap resamples, model refitted in each.
# Output: output/figure5/

source("R/00_setup.R")
df <- load_data()
OUT_DIR <- out_dir("figure5")
ggplot2::theme_set(theme_paper())
set.seed(42)

B_BOOT <- 3000; OUT <- "mort_28d"
GREEN_L <- "#C0E3C0"; GREEN_M <- "#8FCB8F"
L <- df[paired(df), ]; N <- nrow(L)

TG <- tibble::tribble(
  ~v,           ~lab,                  ~cls,
  "rel_bce",    "Coupling trajectory", "Mitochondrial",
  "bce_72h",    "Coupling at 72 h",    "Mitochondrial",
  "ci_0h",      "Complex I",           "Mitochondrial",
  "cii_0h",     "Complex II",          "Mitochondrial",
  "lactate_0h", "Arterial lactate",    "Clinical",
  "d_il10",     "IL-10 trajectory",    "Inflammatory",
  "il10_72h",   "IL-10 at 72 h",       "Inflammatory",
  "tnfa_0h",    "TNF-\u03b1",          "Inflammatory",
  "leptin_0h",  "Leptin",              "Adipokine",
  "lepghre_0h", "Leptin/ghrelin",      "Adipokine",
  "gdf15_0h",   "GDF15",               "Mitokine",
  "fgf21_0h",   "FGF21",               "Mitokine")
CLS <- c(Mitochondrial = TEAL, Clinical = BROWN, Inflammatory = PLUM,
         Adipokine = SAND, Mitokine = ORANG)

coefs <- function(d) {
  m <- tryCatch(suppressWarnings(glm(y ~ x, d, family = binomial)), error = function(e) NULL)
  if (is.null(m) || any(!is.finite(coef(m)))) return(NULL)
  p <- fitted(m)
  if (any(p < 1e-6 | p > 1 - 1e-6)) {
    fi <- tryCatch(logistf::logistf(y ~ x, data = d, pl = FALSE), error = function(e) NULL)
    if (!is.null(fi)) return(list(b = as.numeric(coef(fi)), kind = "Firth", p = NA_real_))
  }
  list(b = as.numeric(coef(m)), kind = "ML", p = summary(m)$coefficients[2, 4])
}
risk <- function(b, x) plogis(b[1] + b[2] * x)
boot_arr <- function(d, shift, B = B_BOOT) vapply(seq_len(B), function(i) {
  db <- d[sample(nrow(d), replace = TRUE), ]
  if (length(unique(db$y)) < 2) return(NA_real_)
  mb <- tryCatch(suppressWarnings(glm(y ~ x, db, family = binomial)), error = function(e) NULL)
  if (is.null(mb) || any(!is.finite(coef(mb)))) return(NA_real_)
  b <- as.numeric(coef(mb))
  mean(risk(b, db$x)) - mean(risk(b, db$x + shift)) }, numeric(1))

# ---- D and E: one SD of improvement ----------------------------------------
ST <- dplyr::bind_rows(lapply(seq_len(nrow(TG)), function(i) {
  d <- tibble::tibble(y = L[[OUT]], x = L[[TG$v[i]]]); d <- d[complete.cases(d), ]
  cf <- coefs(d)
  shift <- ifelse(cf$b[2] > 0, -1, 1) * sd(d$x)
  base <- mean(risk(cf$b, d$x)); arr <- base - mean(risk(cf$b, d$x + shift))
  ic <- quantile(boot_arr(d, shift), c(.025, .975), na.rm = TRUE)
  tibble::tibble(v = TG$v[i], lab = TG$lab[i], cls = TG$cls[i], n = nrow(d), events = sum(d$y),
    baseline = 100 * base, shift = shift, estimator = cf$kind, p_slope = cf$p,
    ARR = 100 * arr, lo = 100 * ic[1], hi = 100 * ic[2], NNT = 1 / arr,
    NNT_lo = ifelse(ic[1] > 0, 1 / ic[2], NA_real_), NNT_hi = ifelse(ic[1] > 0, 1 / ic[1], NA_real_),
    crosses_zero = ic[1] <= 0 & ic[2] >= 0)
})) %>% arrange(desc(ARR))

pD <- ggplot(ST %>% arrange(ARR) %>% mutate(lab = factor(lab, levels = lab),
                                            al = ifelse(crosses_zero, .38, 1)),
             aes(lab, ARR, fill = cls, alpha = al)) +
  geom_hline(yintercept = 0, colour = BLACK, linewidth = .4) +
  geom_col(width = .66, colour = BLACK, linewidth = .28) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .18, linewidth = .36, alpha = 1) +
  scale_fill_manual(values = CLS) + scale_alpha_identity() + coord_flip() +
  labs(x = NULL, y = "Absolute risk reduction (%), one SD of improvement",
       title = "ALL TARGETS, COMMON SCALE") + LEG +
  theme(axis.text.y = element_text(size = 9))
pE <- ggplot(ST %>% filter(!crosses_zero, is.finite(NNT)) %>% arrange(desc(NNT)) %>%
               mutate(lab = factor(lab, levels = lab)), aes(lab, NNT, fill = cls)) +
  geom_col(width = .64, colour = BLACK, linewidth = .28) +
  geom_errorbar(aes(ymin = NNT_lo, ymax = NNT_hi), width = .18, linewidth = .36) +
  scale_fill_manual(values = CLS) + coord_flip() +
  labs(x = NULL, y = "Number needed to treat", title = "TARGETS WITH INTERVAL EXCLUDING ZERO") +
  LEG + theme(axis.text.y = element_text(size = 9))

# ---- A to C: coupling trajectory in natural units --------------------------
dm <- tibble::tibble(y = L[[OUT]], x = L$rel_bce); dm <- dm[complete.cases(dm), ]
cm <- coefs(dm); base_m <- mean(risk(cm$b, dm$x))
GAINS <- c(10, 20, 30)
MT <- dplyr::bind_rows(lapply(GAINS, function(g) {
  arr <- base_m - mean(risk(cm$b, dm$x + g))
  ic <- quantile(boot_arr(dm, g), c(.025, .975), na.rm = TRUE)
  tibble::tibble(gain_pp = g, ARR = 100 * arr, lo = 100 * ic[1], hi = 100 * ic[2],
                 NNT = 1 / arr, NNT_lo = 1 / ic[2], NNT_hi = 1 / ic[1]) })) %>%
  mutate(cen = factor(paste0("+", gain_pp, " pp"), levels = paste0("+", GAINS, " pp")))

VD <- c(Baseline = SAND_L, `+10 pp` = GREEN_L, `+20 pp` = GREEN_M, `+30 pp` = GREEN)
dens <- dplyr::bind_rows(lapply(c(0, GAINS), function(g)
  tibble::tibble(cen = if (g == 0) "Baseline" else paste0("+", g, " pp"),
                 r = 100 * risk(cm$b, dm$x + g)))) %>%
  mutate(cen = factor(cen, levels = names(VD)))
pA <- ggplot(dens, aes(r, fill = cen)) +
  geom_density(alpha = .55, colour = BLACK, linewidth = .25) +
  scale_fill_manual(values = VD) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Predicted 28-day risk (%)", y = "Density", title = "MITOCHONDRIAL TARGET") + LEG +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
bar <- function(ycol, lo, hi, ylab, ttl, top) ggplot(MT, aes(cen, .data[[ycol]], fill = cen)) +
  geom_col(width = .62, colour = BLACK, linewidth = .3) +
  geom_errorbar(aes(ymin = .data[[lo]], ymax = .data[[hi]]), width = .15, linewidth = .36) +
  scale_fill_manual(values = unname(VD[-1])) +
  scale_y_continuous(limits = c(0, top * 1.15), expand = c(0, 0)) +
  labs(x = "Gain in coupling trajectory", y = ylab, title = ttl)
pB <- bar("ARR", "lo", "hi", "Absolute risk reduction (%)", "RISK REDUCTION", max(MT$hi))
pC <- bar("NNT", "NNT_lo", "NNT_hi", "Number needed to treat", "NUMBER NEEDED TO TREAT", max(MT$NNT_hi))

# ---- F: enrichment ----------------------------------------------------------
sel <- dm$x <= 0; n_sel <- sum(sel); GAIN <- 30
r0 <- risk(cm$b, dm$x)
arr_all <- mean(r0) - mean(risk(cm$b, dm$x + GAIN))
arr_sel <- mean(r0) - mean(risk(cm$b, ifelse(sel, dm$x + GAIN, dm$x)))
nnt_all <- 1 / arr_all; nnt_sel <- n_sel / (N * arr_sel)
bs <- t(vapply(seq_len(B_BOOT), function(i) {
  db <- dm[sample(nrow(dm), replace = TRUE), ]
  if (length(unique(db$y)) < 2) return(rep(NA_real_, 4))
  mb <- tryCatch(suppressWarnings(glm(y ~ x, db, family = binomial)), error = function(e) NULL)
  if (is.null(mb) || any(!is.finite(coef(mb)))) return(rep(NA_real_, 4))
  b <- as.numeric(coef(mb)); s <- db$x <= 0; ns <- sum(s)
  if (ns < 3) return(rep(NA_real_, 4))
  a0 <- mean(risk(b, db$x))
  aa <- a0 - mean(risk(b, db$x + GAIN)); as <- a0 - mean(risk(b, ifelse(s, db$x + GAIN, db$x)))
  c(aa, as, 1 / aa, ns / (nrow(db) * as)) }, numeric(4)))
q <- function(j) quantile(bs[, j], c(.025, .975), na.rm = TRUE)
EN <- tibble::tibble(strategy = factor(c("TREAT ALL", "TREAT SELECTED"), levels = c("TREAT ALL", "TREAT SELECTED")),
  exposed = c(N, n_sel), ARR = 100 * c(arr_all, arr_sel), NNT = c(nnt_all, nnt_sel),
  lo = c(q(3)[1], q(4)[1]), hi = c(q(3)[2], q(4)[2])) %>%
  mutate(lab = factor(sprintf("%s\n(%d exposed)", c("Treat all", "Treat selected"), exposed),
                      levels = sprintf("%s\n(%d exposed)", c("Treat all", "Treat selected"), c(N, n_sel))))
pF <- ggplot(EN, aes(lab, NNT, fill = strategy)) +
  geom_col(width = .5, colour = BLACK, linewidth = .3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .14, linewidth = .36) +
  scale_fill_manual(values = c(`TREAT ALL` = SAND, `TREAT SELECTED` = GREEN)) +
  scale_y_continuous(limits = c(0, max(EN$hi, na.rm = TRUE) * 1.12), expand = c(0, 0)) +
  labs(x = NULL, y = "Number needed to treat, per patient treated",
       title = "ENRICHMENT BY DIRECT MEASUREMENT")

# ---- save -------------------------------------------------------------------
fig5 <- (pA | pB | pC) / (pD | pE | pF) + patchwork::plot_annotation(tag_levels = "A") & TAG
save_fig(fig5, OUT_DIR, "Figure5", 13.6, 8.6)
write_out(MT %>% select(-cen), OUT_DIR, "panels_A_C_coupling_gains.csv")
write_out(ST, OUT_DIR, "panels_D_E_targets.csv")
write_out(EN %>% select(-lab) %>% mutate(benefit_retained_pct = c(100, 100 * arr_sel / arr_all)),
          OUT_DIR, "panel_F_enrichment.csv")
writeLines(c(sprintf("Patients: %d; 28-day deaths: %d", N, sum(L[[OUT]])),
             sprintf("Baseline predicted risk: %.1f%% (%s estimator)", 100 * base_m, cm$kind),
             sprintf("Targets with interval excluding zero: %d of %d", sum(!ST$crosses_zero), nrow(ST))),
           file.path(OUT_DIR, "figure5_summary.txt"))
