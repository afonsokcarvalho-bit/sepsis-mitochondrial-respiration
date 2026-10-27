# Table 3. Nested model comparisons in the 37 patients with paired
# respirometry: penalised likelihood ratio test (logistftest, added
# coefficients constrained to zero within the full Firth model, standardised
# predictors) and cross-validated gain in AUC (five folds, 20 repetitions,
# Firth inside each training fold). Seven comparisons x four horizons = 28
# tests; Bonferroni threshold 0.05 / 28.
# Output: output/table3/Table3.csv and table3_long.csv

source("R/00_setup.R")
df <- load_data()
OUT <- out_dir("table3")
set.seed(42)
K_CV <- 5; R_CV <- 20

ORD <- by_events(df)
L <- df[paired(df), ]
RAW <- c(z_r = "r_0h", z_ci = "ci_0h", z_cii = "cii_0h", z_bce = "bce_0h",
         z_traj = "rel_bce", z_sofa1 = "sofa_0h", z_sofa3 = "sofa_72h")
for (k in names(RAW)) L[[k]] <- as.numeric(scale(L[[RAW[[k]]]]))

auc_of <- function(y, p) if (length(unique(y)) < 2) NA_real_ else
  as.numeric(pROC::auc(pROC::roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)))
ffit <- function(f, d) tryCatch(logistf::logistf(f, data = d, pl = FALSE,
  control = logistf::logistf.control(maxit = 500, maxstep = 5)), error = function(e) NULL)
pred_f <- function(m, f, d) { X <- model.matrix(f, d); b <- coef(m)[colnames(X)]
  as.numeric(plogis(as.numeric(X %*% b))) }
cv_firth <- function(vs, d, y) {
  f <- as.formula(paste(y, "~", paste(vs, collapse = "+")))
  n <- nrow(d); out <- numeric(0)
  for (r in seq_len(R_CV)) {
    fold <- sample(rep(seq_len(K_CV), length.out = n))
    p <- rep(NA_real_, n); bad <- FALSE
    for (k in seq_len(K_CV)) {
      tr <- fold != k
      if (length(unique(d[[y]][tr])) < 2) { bad <- TRUE; break }
      m <- ffit(f, d[tr, , drop = FALSE]); if (is.null(m)) { bad <- TRUE; break }
      pk <- tryCatch(pred_f(m, f, d[!tr, , drop = FALSE]), error = function(e) NULL)
      if (is.null(pk) || any(!is.finite(pk))) { bad <- TRUE; break }
      p[!tr] <- pk
    }
    if (bad || any(is.na(p))) next
    out <- c(out, auc_of(d[[y]], p))
  }
  if (length(out)) mean(out) else NA_real_
}
p_penalised <- function(full, added, d, y) {
  m <- ffit(as.formula(paste(y, "~", paste(full, collapse = "+"))), d)
  if (is.null(m)) return(NA_real_)
  t <- tryCatch(logistf::logistftest(m, test = as.formula(
         paste("~", paste(added, collapse = "+"), "- 1"))), error = function(e) NULL)
  if (is.null(t)) NA_real_ else as.numeric(t$prob)
}

NEST <- list(
  list(f = c("z_ci", "z_traj"),    b = "z_ci",    q = "Trajectory added to Complex I"),
  list(f = c("z_ci", "z_traj"),    b = "z_traj",  q = "Complex I added to trajectory"),
  list(f = c("z_r", "z_ci", "z_cii", "z_bce", "z_traj"),
       b = c("z_r", "z_ci", "z_cii", "z_bce"),    q = "Trajectory added to admission profile"),
  list(f = c("z_traj", "z_sofa1"), b = "z_sofa1", q = "Trajectory added to SOFA, 0 h"),
  list(f = c("z_traj", "z_sofa1"), b = "z_traj",  q = "SOFA, 0 h, added to trajectory"),
  list(f = c("z_traj", "z_sofa3"), b = "z_sofa3", q = "Trajectory added to SOFA, 72 h"),
  list(f = c("z_ci", "z_sofa1"),   b = "z_sofa1", q = "Complex I added to SOFA, 0 h"))

RES <- dplyr::bind_rows(lapply(ORD, function(y) dplyr::bind_rows(lapply(NEST, function(z) {
  d <- L[, c(y, z$f)]; d <- d[complete.cases(d), ]
  cb <- cv_firth(z$b, d, y); cf <- cv_firth(z$f, d, y)
  tibble::tibble(horizon = nice(y), comparison = z$q, n = nrow(d), events = sum(d[[y]]),
                 P = p_penalised(z$f, setdiff(z$f, z$b), d, y),
                 auc_smaller = cb, auc_larger = cf, gain = cf - cb) }))))
write_out(RES, OUT, "table3_long.csv")

cell <- function(p, g) sprintf("%s (%+.3f)", ifelse(p < 1e-4, "< 0.0001",
  ifelse(p >= 0.01, formatC(signif(p, 2), format = "fg", digits = 2, flag = "#"),
         formatC(signif(p, 2), format = "fg", digits = 2))), g)
TAB3 <- RES %>% mutate(cell = cell(P, gain)) %>% select(comparison, horizon, cell) %>%
  pivot_wider(names_from = horizon, values_from = cell) %>%
  dplyr::slice(match(vapply(NEST, `[[`, "", "q"), comparison)) %>%
  select(Comparison = comparison, dplyr::all_of(HORIZONS))
write_out(TAB3, OUT, "Table3.csv")
message("Bonferroni threshold: P < ", signif(0.05 / nrow(RES), 2))
