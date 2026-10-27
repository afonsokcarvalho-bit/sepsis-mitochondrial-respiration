# Figure 4 and fig. S1.
#   A  Principal components of the four admission respirometry parameters
#      (log1p, standardised; n = 40); circulating and clinical markers projected
#      as supplementary variables by Spearman correlation with the component
#      scores (Benjamini-Hochberg across markers, per component).
#   B, C  ROC curves for 28-day mortality, from pooled out-of-fold predictions.
#   D, E  Cross-validated AUC with 95% bootstrap intervals, 28 days and 6 months.
#   F  Mortality by baseline coupling (median split) and direction of change,
#      Wilson intervals, Fisher exact tests.
#   fig. S1  Coupling trajectories, with the interquartile range of six
#      healthy controls assayed with the same protocol (ref. 41).
# Models in B to E are fitted and evaluated on the 37 patients with paired
# respirometry: five folds, 20 repetitions, Firth logistic regression inside
# each training fold. The interval is the bootstrap percentile interval
# (1,500 resamples) of the mean per-repetition AUC.
# Output: output/figure4/

source("R/00_setup.R")
df <- load_data()
OUT_DIR <- out_dir("figure4")
ggplot2::theme_set(theme_paper())
set.seed(42)

K_CV <- 5; R_CV <- 20; OUT <- "mort_28d"
LM_KEEP <- paired(df)

# ---- cross-validation engine ------------------------------------------------
auc_of <- function(y, p) if (length(unique(y)) < 2 || all(!is.finite(p))) NA_real_ else
  as.numeric(pROC::auc(pROC::roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)))

fit_firth <- function(vs, d, y) {
  f <- as.formula(paste(y, "~", paste(vs, collapse = "+")))
  m <- tryCatch(logistf::logistf(f, data = d, pl = FALSE,
         control = logistf::logistf.control(maxit = 500, maxstep = 5)),
       error = function(e) NULL)
  if (is.null(m)) return(NULL)
  list(m = m, f = f)
}
pred_firth <- function(o, d) {
  X <- model.matrix(o$f, d); b <- coef(o$m)[colnames(X)]
  as.numeric(plogis(as.numeric(X %*% b)))
}
fauc <- function(yy, pp) { n1 <- sum(yy == 1); n0 <- sum(yy == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(pp)[yy == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }

cv_pack <- function(vs, y = OUT, K = K_CV, R = R_CV, keep = LM_KEEP) {
  d <- df[keep, c(y, vs)]; d <- d[complete.cases(d), ]
  if (nrow(d) < 18 || length(unique(d[[y]])) < 2) return(NULL)
  n <- nrow(d); P <- matrix(NA_real_, n, R); A <- rep(NA_real_, R)
  for (r in seq_len(R)) {
    fold <- sample(rep(seq_len(K), length.out = n))
    p <- rep(NA_real_, n); bad <- FALSE
    for (k in seq_len(K)) {
      tr <- fold != k
      if (length(unique(d[[y]][tr])) < 2) { bad <- TRUE; break }
      o <- fit_firth(vs, d[tr, ], y)
      if (is.null(o)) { bad <- TRUE; break }
      pk <- tryCatch(pred_firth(o, d[!tr, , drop = FALSE]), error = function(e) NULL)
      if (is.null(pk) || any(!is.finite(pk))) { bad <- TRUE; break }
      p[!tr] <- pk
    }
    if (bad || any(is.na(p))) next
    P[, r] <- p; A[r] <- auc_of(d[[y]], p)
  }
  ok <- !is.na(A)
  if (!any(ok)) return(NULL)
  pm <- rowMeans(P[, ok, drop = FALSE], na.rm = TRUE)
  PK <- P[, ok, drop = FALSE]
  bs <- replicate(1500, {
    i <- sample(n, replace = TRUE); yy <- d[[y]][i]
    if (length(unique(yy)) < 2) NA_real_ else
      mean(apply(PK[i, , drop = FALSE], 2, function(pp) fauc(yy, pp))) })
  ic <- quantile(bs, c(.025, .975), na.rm = TRUE)
  list(auc = mean(A[ok]), sd = sd(A[ok]), lo = min(ic[1], 1), hi = min(ic[2], 1),
       reps = sum(ok), n = n, events = sum(d[[y]]), y = d[[y]], p = pm,
       roc = pROC::roc(d[[y]], pm, levels = c(0, 1), direction = "<", quiet = TRUE))
}

# ---- models -----------------------------------------------------------------
RESP <- list(
  list(id = "R1", lab = "Admission respiratory profile",  v = c("r_0h","ci_0h","cii_0h","bce_0h")),
  list(id = "T1", lab = "Coupling trajectory",            v = "rel_bce"),
  list(id = "T2", lab = "Respiratory profile trajectory", v = c("rel_r","rel_ci","rel_cii","rel_bce")),
  list(id = "R5", lab = "Respiratory profile at 72 h",    v = c("r_72h","ci_72h","cii_72h","bce_72h")),
  list(id = "R2", lab = "Complex I alone",                v = "ci_0h"),
  list(id = "R3", lab = "Coupling at admission",          v = "bce_0h"))
CIRC <- list(
  list(id = "C1", lab = "Organ-failure score",            v = "sofa_0h"),
  list(id = "C3", lab = "Lactate + score",                v = c("lactate_0h","sofa_0h")),
  list(id = "C4", lab = "Organ-failure score, day 3",     v = "sofa_72h"),
  list(id = "I2", lab = "IL-10 trajectory",               v = "d_il10"),
  list(id = "I1", lab = "IL-10 72 h + TNF-\u03b1",         v = c("il10_72h","tnfa_0h")),
  list(id = "M1", lab = "Stress mitokines",               v = c("gdf15_0h","fgf21_0h")),
  list(id = "M2", lab = "Mitokines + adipokines",         v = c("gdf15_0h","fgf21_0h","leptin_0h","ghrelin_0h")))
COMB <- list(                                   # X8 and X9 are exploratory
  list(id = "X6", lab = "Coupling trajectory + score",    v = c("rel_bce","sofa_0h")),
  list(id = "X8", lab = "Coupling + IL-10 trajectories",  v = c("rel_bce","d_il10")),
  list(id = "X9", lab = "Both trajectories + score",      v = c("rel_bce","d_il10","sofa_0h")))
ALL <- c(RESP, CIRC, COMB)
ids <- function(L) vapply(L, function(z) z$id, "")

run_models <- function(y) {
  fits <- list()
  for (z in ALL) {
    s <- cv_pack(z$v, y = y)
    if (!is.null(s)) { s$lab <- z$lab; s$id <- z$id; fits[[z$id]] <- s }
  }
  fits
}
auc_table <- function(fits) dplyr::bind_rows(lapply(names(fits), function(k) {
  s <- fits[[k]]
  cls <- if (k %in% ids(RESP)) "Respiratory" else if (k %in% ids(COMB)) "Combined"
         else "Circulating and clinical"
  tibble::tibble(id = k, lab = s$lab, cls = cls, n = s$n, events = s$events,
                 auc = s$auc, sd = s$sd, lo = s$lo, hi = s$hi, reps = s$reps) }))

# 28-day models are fitted first, then the same models at 6 months; this
# order reproduces the published values.
FIT  <- run_models(OUT)
AUCS <- auc_table(FIT); AUC_REF <- AUCS$auc[AUCS$id == "C1"]
CLS <- c(Respiratory = TEAL, `Circulating and clinical` = PLUM, Combined = GREEN)

# ---- panel A ----------------------------------------------------------------
ACT  <- c("r_0h","ci_0h","cii_0h","bce_0h")
ALAB <- c(r_0h = "LEAK", ci_0h = "CI", cii_0h = "CII", bce_0h = "Coupling")
dpc <- df[, c(".row", ACT, OUT)]; dpc <- dpc[complete.cases(dpc), ]
M <- scale(as.data.frame(lapply(dpc[, ACT], function(x) log1p(pmax(x, 0)))))
pca <- prcomp(M); imp <- summary(pca)$importance
sc <- as.data.frame(pca$x[, 1:2]); names(sc) <- c("PC1","PC2")
sc$.row <- dpc$.row
sc$grp <- factor(dpc[[OUT]], labels = c("SURVIVOR","NON-SURVIVOR"))
smax <- max(abs(c(sc$PC1, sc$PC2)))
act <- as.data.frame(pca$rotation[, 1:2]) %>% tibble::rownames_to_column("v") %>%
  mutate(lab = unname(ALAB[v]), x = PC1 * smax * .92, y = PC2 * smax * .92)

SUPP <- c(lactate_0h = "Lactate", sofa_0h = "SOFA", il10_72h = "IL-10 72h",
          d_il10 = "IL-10 change", tnfa_0h = "TNF-\u03b1", leptin_0h = "Leptin",
          lepghre_0h = "Leptin/ghrelin", gdf15_0h = "GDF15", fgf21_0h = "FGF21")
proj <- dplyr::bind_rows(lapply(names(SUPP), function(v) {
  d <- dplyr::inner_join(sc, df[, c(".row", v)], by = ".row")
  names(d)[ncol(d)] <- "x"; ok <- !is.na(d$x)
  if (sum(ok) < 15) return(NULL)
  c1 <- suppressWarnings(cor.test(d$PC1[ok], d$x[ok], method = "spearman"))
  c2 <- suppressWarnings(cor.test(d$PC2[ok], d$x[ok], method = "spearman"))
  tibble::tibble(lab = unname(SUPP[v]), n = sum(ok),
                 rho1 = as.numeric(c1$estimate), p1 = c1$p.value,
                 rho2 = as.numeric(c2$estimate), p2 = c2$p.value) })) %>%
  mutate(fdr1 = p.adjust(p1, "BH"), fdr2 = p.adjust(p2, "BH"),
         norm = sqrt(rho1^2 + rho2^2), x = rho1 * smax * .92, y = rho2 * smax * .92,
         sig = fdr1 < .05 | fdr2 < .05)

pA <- ggplot() +
  geom_hline(yintercept = 0, colour = GREY, linewidth = .3, linetype = "dotted") +
  geom_vline(xintercept = 0, colour = GREY, linewidth = .3, linetype = "dotted") +
  stat_ellipse(data = sc, aes(PC1, PC2, colour = grp), level = .68, linewidth = .45,
               linetype = "dashed", show.legend = FALSE) +
  geom_point(data = sc, aes(PC1, PC2, fill = grp), shape = 21, size = 2.5,
             colour = BLACK, stroke = .3) +
  geom_segment(data = proj, aes(x = 0, y = 0, xend = x, yend = y,
               linetype = ifelse(sig, "s", "n"), alpha = ifelse(sig, 1, .45)),
               colour = PLUM, linewidth = .5, show.legend = FALSE,
               arrow = arrow(length = unit(.14, "cm"), type = "closed")) +
  ggrepel::geom_text_repel(data = proj, aes(x, y, label = lab, alpha = ifelse(sig, 1, .5),
               fontface = ifelse(sig, "bold", "plain")), colour = PLUM, family = FONT,
               size = 3, segment.size = .18, min.segment.length = .3, max.overlaps = 40,
               show.legend = FALSE) +
  scale_alpha_identity() +
  geom_segment(data = act, aes(x = 0, y = 0, xend = x, yend = y), colour = TEAL,
               linewidth = .8, arrow = arrow(length = unit(.17, "cm"), type = "closed")) +
  ggrepel::geom_text_repel(data = act, aes(x, y, label = lab), colour = TEAL,
               family = FONT, size = 3.4, fontface = "bold", segment.size = .18,
               max.overlaps = 40) +
  scale_linetype_manual(values = c(s = "solid", n = "22")) +
  scale_fill_manual(values = c(SURVIVOR = GREEN, `NON-SURVIVOR` = ORANG)) +
  scale_colour_manual(values = c(SURVIVOR = GREEN, `NON-SURVIVOR` = ORANG)) +
  labs(x = sprintf("PC1 (%.1f%%)", 100 * imp[2, 1]), y = sprintf("PC2 (%.1f%%)", 100 * imp[2, 2]),
       title = "BIOENERGETIC SPACE, MARKERS PROJECTED") + LEG

# ---- panels B and C ---------------------------------------------------------
roc_xy <- function(s) tibble::tibble(x = 1 - rev(s$roc$specificities), y = rev(s$roc$sensitivities))
ref_d <- roc_xy(FIT[["C1"]])
roc_panel <- function(sel, pal, ttl) {
  d <- dplyr::bind_rows(lapply(sel, function(k) { s <- FIT[[k]]; if (is.null(s)) return(NULL)
    roc_xy(s) %>% mutate(lab = sprintf("%s  (%.2f)", s$lab, s$auc), ord = s$auc) }))
  d$lab <- factor(d$lab, levels = unique(d$lab[order(-d$ord)]))
  ggplot(d, aes(x, y, colour = lab)) +
    geom_abline(linetype = "dotted", colour = GREY, linewidth = .32) +
    geom_line(data = ref_d, aes(x, y), colour = REF, linewidth = .75,
              linetype = "longdash", inherit.aes = FALSE) +
    geom_line(linewidth = .68) +
    scale_colour_manual(values = pal[seq_len(dplyr::n_distinct(d$lab))]) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
    annotate("text", x = .98, y = .045, label = sprintf("SOFA  (%.2f)", AUC_REF),
             hjust = 1, vjust = 0, family = FONT, size = 3, fontface = "bold", colour = REF) +
    labs(x = "1 - specificity", y = "Sensitivity", title = ttl) +
    theme(legend.position = c(.66, .27), legend.title = element_blank(),
          legend.text = element_text(size = 8.2, family = FONT),
          legend.key.size = unit(.32, "cm"), legend.key.width = unit(.68, "cm"),
          legend.background = element_rect(fill = "white", colour = NA))
}
pB <- roc_panel(c(ids(RESP), "X6", "X9"),
                c("#0F4A52", TEAL, "#5FB0AE", GREEN, "#2E6B2E", SAND, "#8FCB8F", "#1F6B6B"),
                "RESPIRATORY MODELS")
pC <- roc_panel(setdiff(ids(CIRC), "C1"),
                c("#4A3A57", PLUM, "#B79FC2", SAND, ORANG, "#B5651D"),
                "CIRCULATING AND CLINICAL MODELS")

# ---- same models, 6-month mortality ----------------------------------------
FIT6  <- run_models("mort_6m")
AUCS6 <- auc_table(FIT6); AUC_REF6 <- AUCS6$auc[AUCS6$id == "C1"]

# ---- panels D and E ---------------------------------------------------------
rank_panel <- function(A, ref, ttl, legend = TRUE) {
  ad <- A %>% arrange(auc) %>% mutate(lab = factor(lab, levels = lab), truncated = hi >= .9995)
  p <- ggplot(ad, aes(lab, auc, colour = cls)) +
    geom_hline(yintercept = .5, linetype = "dotted", colour = GREY, linewidth = .35) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = ref, ymax = 1, fill = REF, alpha = .05) +
    geom_hline(yintercept = ref, linetype = "longdash", colour = REF, linewidth = .55) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = .15, linewidth = .42) +
    geom_point(aes(fill = cls), shape = 21, size = 3, colour = BLACK, stroke = .32) +
    geom_point(data = ad[ad$truncated, ], aes(y = 1), shape = 62, size = 2.1,
               colour = BLACK, show.legend = FALSE) +
    scale_colour_manual(values = CLS) + scale_fill_manual(values = CLS) +
    scale_y_continuous(limits = c(.30, 1.02), breaks = seq(.3, 1, .1), oob = scales::squish) +
    coord_flip() +
    labs(x = NULL, y = "Cross-validated area (AUC)", title = ttl) + LEG +
    theme(axis.text.y = element_text(size = 8.6))
  if (!legend) p <- p + theme(legend.position = "none")
  p
}
pD <- rank_panel(AUCS,  AUC_REF,  "DISCRIMINATION, 28 DAYS")
pE <- rank_panel(AUCS6, AUC_REF6, "DISCRIMINATION, 6 MONTHS", legend = FALSE)

# ---- fig. S1 ----------------------------------------------------------------
dtr <- df %>% filter(!is.na(d_bce))
med <- dtr %>% group_by(mort_28d) %>%
  summarise(b0 = median(bce_0h), b72 = median(bce_72h), .groups = "drop") %>%
  mutate(g = factor(mort_28d, labels = c("SURVIVOR","NON-SURVIVOR")))
HEALTHY_LO <- 0.40; HEALTHY_HI <- 0.50     # IQR, six healthy controls (ref. 41)
pS1 <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = HEALTHY_LO, ymax = HEALTHY_HI,
           fill = SAND_L, alpha = .55) +
  annotate("text", x = 36, y = HEALTHY_HI + .035, label = "HEALTHY CONTROLS (n = 6), IQR",
           family = FONT, size = 2.9, colour = "#8A6F1A", fontface = "italic") +
  geom_segment(data = dtr, aes(x = 0, xend = 72, y = bce_0h, yend = bce_72h, group = .row),
               colour = GREY, alpha = .32, linewidth = .3) +
  geom_segment(data = med, aes(x = 0, xend = 72, y = b0, yend = b72, colour = g), linewidth = 1.5) +
  geom_point(data = med, aes(0, b0, fill = g), shape = 21, size = 3, colour = BLACK,
             stroke = .38, show.legend = FALSE) +
  geom_point(data = med, aes(72, b72, fill = g), shape = 21, size = 3, colour = BLACK,
             stroke = .38, show.legend = FALSE) +
  scale_colour_manual(values = c(SURVIVOR = GREEN, `NON-SURVIVOR` = ORANG)) +
  scale_fill_manual(values = c(SURVIVOR = GREEN, `NON-SURVIVOR` = ORANG), guide = "none") +
  scale_x_continuous(breaks = c(0, 72), labels = c("0 h", "72 h"), limits = c(-7, 79),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, .95), expand = c(0, 0)) +
  labs(x = NULL, y = "(P-Ln)/P", title = "COUPLING TRAJECTORY, 28-DAY OUTCOME") + LEG +
  theme(axis.text.x = element_text(face = "bold", size = 11))

# ---- panel F ----------------------------------------------------------------
cut0 <- median(dtr$bce_0h, na.rm = TRUE)
dtr <- dtr %>% mutate(level = ifelse(bce_0h < cut0, "Low", "High"),
                      direction = ifelse(d_bce > 0, "recovered", "lost"))
g4 <- dtr %>% mutate(gr = paste0(level, "\n", direction)) %>% group_by(gr) %>%
  summarise(n = dplyr::n(), o28 = sum(mort_28d), o6 = sum(mort_6m), .groups = "drop") %>%
  pivot_longer(c(o28, o6), names_to = "hz", values_to = "ob") %>%
  mutate(hz = factor(ifelse(hz == "o28", "28 DAYS", "6 MONTHS"), levels = c("28 DAYS", "6 MONTHS")),
         pct = 100 * ob / n) %>%
  rowwise() %>% mutate(lo = wilson(ob, n)[1], hi = wilson(ob, n)[2]) %>% ungroup() %>%
  mutate(gr = factor(gr, levels = c("High\nlost", "Low\nlost", "High\nrecovered", "Low\nrecovered")))
pF <- ggplot(g4, aes(gr, pct, fill = hz)) +
  geom_col(position = position_dodge(.78), width = .66, colour = BLACK, linewidth = .3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .15, position = position_dodge(.78), linewidth = .35) +
  scale_fill_manual(values = c(`28 DAYS` = GREEN, `6 MONTHS` = ORANG)) +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = NULL, y = "Mortality (%)", title = "BASELINE AND TRAJECTORY") + LEG +
  theme(axis.text.x = element_text(size = 9))

fisher <- function(g, y) fisher.test(table(g, y))$p.value
FT <- tibble::tibble(
  comparison = c("four groups", "direction (lost vs recovered)", "level (high vs low)"),
  P_28d = c(fisher(paste(dtr$level, dtr$direction), dtr$mort_28d),
            fisher(dtr$d_bce <= 0, dtr$mort_28d), fisher(dtr$bce_0h >= cut0, dtr$mort_28d)),
  P_6m  = c(fisher(paste(dtr$level, dtr$direction), dtr$mort_6m),
            fisher(dtr$d_bce <= 0, dtr$mort_6m),  fisher(dtr$bce_0h >= cut0, dtr$mort_6m)))
POOL <- dtr %>% group_by(direction) %>%
  summarise(n = dplyr::n(), mort_28d = 100 * mean(mort_28d), mort_6m = 100 * mean(mort_6m),
            .groups = "drop")

# ---- save -------------------------------------------------------------------
fig4 <- (pA | pB | pC) / (pD | pE | pF) + patchwork::plot_annotation(tag_levels = "A") & TAG
save_fig(fig4, OUT_DIR, "Figure4", 13.6, 8.6)
save_fig(pS1, OUT_DIR, "FigureS1", 5.2, 4.6)

pc <- as.data.frame(pca$rotation[, 1:2]) %>% tibble::rownames_to_column("parameter") %>%
  mutate(variance_PC1 = imp[2, 1], variance_PC2 = imp[2, 2])
write_out(pc, OUT_DIR, "panelA_pca_loadings.csv")
write_out(proj %>% select(marker = lab, n, rho1, fdr1, rho2, fdr2), OUT_DIR, "panelA_projection.csv")
write_out(AUCS,  OUT_DIR, "panelD_auc_28d.csv")
write_out(AUCS6, OUT_DIR, "panelE_auc_6m.csv")
write_out(g4 %>% mutate(gr = gsub("\n", " ", gr)), OUT_DIR, "panelF_groups.csv")
write_out(POOL, OUT_DIR, "panelF_by_direction.csv")
write_out(FT, OUT_DIR, "panelF_fisher.csv")
write_out(med %>% select(group = g, median_0h = b0, median_72h = b72), OUT_DIR, "figS1_medians.csv")
writeLines(c(
  sprintf("Paired respirometry: %d patients", sum(LM_KEEP)),
  sprintf("SOFA reference, 28 days: AUC %.3f; models above it: %d of %d",
          AUC_REF, sum(AUCS$auc > AUC_REF), nrow(AUCS) - 1),
  sprintf("SOFA reference, 6 months: AUC %.3f; models above it: %d of %d",
          AUC_REF6, sum(AUCS6$auc > AUC_REF6), nrow(AUCS6) - 1),
  sprintf("Median baseline coupling (panel F cut): %.4f", cut0)),
  file.path(OUT_DIR, "figure4_summary.txt"))
