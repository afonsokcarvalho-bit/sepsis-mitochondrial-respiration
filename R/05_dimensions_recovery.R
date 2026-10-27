# Respiratory capacity and coupling as two dimensions, and recovery of coupling.
#   1. Spearman correlation of Complex I with coupling efficiency at admission
#      (bootstrap 95% interval, 2,000 resamples) and with Complex II (n = 40).
#   2. Recovery of coupling (any increase in BCE from 0 h to 72 h, n = 37):
#      AUC of admission Complex I (DeLong interval), Spearman correlation with
#      the change and the relative change in coupling, and Firth logistic
#      regression of recovery on standardised Complex I and SOFA.
#   3. Firth coefficient of standardised Complex I for each mortality horizon,
#      alone and adjusted for the standardised coupling trajectory, with a
#      bootstrap interval (1,000 resamples) for the difference. Descriptive;
#      not a mediation analysis.
# Output: output/dimensions_recovery/

source("R/00_setup.R")
df <- load_data()
OUT <- out_dir("dimensions_recovery")
set.seed(42)

L <- df[paired(df), ]
firth_fit <- function(f, d) tryCatch(logistf::logistf(f, data = d, pl = FALSE,
  control = logistf::logistf.control(maxit = 500, maxstep = 5)), error = function(e) NULL)

rho_ci <- function(x, y, B = 2000) {
  ok <- complete.cases(x, y); x <- x[ok]; y <- y[ok]
  r  <- suppressWarnings(cor(x, y, method = "spearman"))
  p  <- suppressWarnings(cor.test(x, y, method = "spearman")$p.value)
  bs <- replicate(B, { i <- sample(length(x), replace = TRUE)
    suppressWarnings(cor(x[i], y[i], method = "spearman")) })
  q <- quantile(bs, c(.025, .975), na.rm = TRUE)
  tibble::tibble(n = length(x), rho = r, lo = unname(q[1]), hi = unname(q[2]), P = p)
}
rho <- function(x, y) { ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
  tibble::tibble(n = sum(complete.cases(x, y)), rho = unname(ct$estimate), P = ct$p.value) }

# ---- 1. two dimensions ------------------------------------------------------
DIM <- dplyr::bind_rows(
  rho_ci(df$ci_0h, df$bce_0h) %>% mutate(pair = "Complex I vs coupling", .before = 1),
  rho(df$ci_0h, df$cii_0h)    %>% mutate(pair = "Complex I vs Complex II", .before = 1))
write_out(DIM, OUT, "capacity_vs_coupling.csv")

# ---- 2. recovery of coupling ------------------------------------------------
L$recovered <- as.integer(L$d_bce > 0)
r_ci <- pROC::roc(L$recovered, L$ci_0h, levels = c(0, 1), direction = "<", quiet = TRUE)
ci_auc <- as.numeric(pROC::ci.auc(r_ci))
REC_AUC <- tibble::tibble(predictor = "Complex I, 0 h", n = length(r_ci$response),
  recovered = sum(L$recovered), AUC = ci_auc[2], lo = ci_auc[1], hi = ci_auc[3])
REC_RHO <- dplyr::bind_rows(
  rho(L$ci_0h, L$d_bce)   %>% mutate(pair = "Complex I, 0 h vs change in coupling", .before = 1),
  rho(L$ci_0h, L$rel_bce) %>% mutate(pair = "Complex I, 0 h vs coupling trajectory", .before = 1))
dr <- L[, c("recovered", "ci_0h", "sofa_0h")]; dr <- dr[complete.cases(dr), ]
dr$zc <- as.numeric(scale(dr$ci_0h)); dr$zs <- as.numeric(scale(dr$sofa_0h))
mr <- firth_fit(recovered ~ zc + zs, dr)
REC_OR <- tibble::tibble(term = c("Complex I, per SD", "SOFA, per SD"),
  OR = exp(coef(mr)[c("zc", "zs")]), lo = exp(mr$ci.lower[c("zc", "zs")]),
  hi = exp(mr$ci.upper[c("zc", "zs")]), P = mr$prob[c("zc", "zs")])
write_out(REC_AUC, OUT, "recovery_auc.csv")
write_out(REC_RHO, OUT, "recovery_correlations.csv")
write_out(REC_OR,  OUT, "recovery_firth.csv")

# ---- 3. Complex I alone and adjusted for the coupling trajectory ------------
set.seed(42)
B_BOOT <- 1000
AT <- dplyr::bind_rows(lapply(OUTCOMES, function(y) {
  d <- L[, c(y, "ci_0h", "rel_bce")]; d <- d[complete.cases(d), ]
  d$zc <- as.numeric(scale(d$ci_0h)); d$zt <- as.numeric(scale(d$rel_bce))
  fa <- as.formula(paste(y, "~ zc")); fb <- as.formula(paste(y, "~ zc + zt"))
  ma <- firth_fit(fa, d); mb <- firth_fit(fb, d)
  bs <- t(replicate(B_BOOT, {
    db <- d[sample(nrow(d), replace = TRUE), ]
    if (length(unique(db[[y]])) < 2) return(c(NA_real_, NA_real_))
    m1 <- firth_fit(fa, db); m2 <- firth_fit(fb, db)
    if (is.null(m1) || is.null(m2)) return(c(NA_real_, NA_real_))
    c(unname(coef(m1)["zc"]), unname(coef(m2)["zc"])) }))
  dif <- bs[, 1] - bs[, 2]
  tibble::tibble(horizon = nice(y), n = nrow(d), events = sum(d[[y]]),
    complex_I_alone = unname(coef(ma)["zc"]), complex_I_adjusted = unname(coef(mb)["zc"]),
    trajectory_adjusted = unname(coef(mb)["zt"]),
    difference_lo = quantile(dif, .025, na.rm = TRUE),
    difference_hi = quantile(dif, .975, na.rm = TRUE))
}))
write_out(AT, OUT, "complex_I_adjusted_for_trajectory.csv")
