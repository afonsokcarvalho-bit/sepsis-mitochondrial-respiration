# Table 2 and tables S1 to S3.
#   table S1  Complex I-linked respiration and coupling efficiency, alone and
#             combined, at 0 h, at 72 h and as trajectories (n = 37):
#               A  regression to the mean, flux changes and ANCOVA
#               B  cross-validated AUC
#               C  comparison with the SOFA score (DeLong on out-of-fold
#                  predictions; penalised likelihood ratio test; gain in AUC)
#               D  one respiratory parameter added to the other
#   table S2  variables by recovery of coupling efficiency (n = 37)
#   table S3  admission models in the 37 patients with paired respirometry
#             and in all 40 patients
#   Table 2   rows of tables S1B and S3 shown in the main text
# Cross-validation: five folds, 20 repetitions, Firth logistic regression with
# standardised predictors inside each training fold. Horizons are processed
# in increasing order of events and the sections in the order below; changing
# either changes the random folds and therefore the published values.
# Output: output/supplementary_tables/ (Supplementary_Tables.docx and .txt,
#         Table2.csv)

source("R/00_setup.R")
df <- load_data()
OUT <- out_dir("supplementary_tables")
set.seed(42)
K_CV <- 5; R_CV <- 20

ORD <- by_events(df)
PAIRED <- !is.na(df$bce_0h) & !is.na(df$bce_72h) & !is.na(df$ci_0h) & !is.na(df$ci_72h)

f3 <- function(x) fmt(x, 3)
fg <- function(x) ifelse(is.na(x), "NA", sprintf("%+.3f", x))

# ---- models (order fixes the random folds) ----------------------------------
MOD <- tibble::tribble(
  ~id,     ~vars,
  "SOFA1", list("sofa_0h"),
  "SOFA3", list("sofa_72h"),
  "CI1",   list("ci_0h"),
  "BCE1",  list("bce_0h"),
  "D1",    list(c("ci_0h", "bce_0h")),
  "CI3",   list("ci_72h"),
  "BCE3",  list("bce_72h"),
  "D3",    list(c("ci_72h", "bce_72h")),
  "TCI",   list("rel_ci"),
  "TCIa",  list("d_ci"),
  "TBCE",  list("rel_bce"),
  "TT",    list(c("rel_ci", "rel_bce")),
  "D13",   list(c("ci_0h", "bce_0h", "ci_72h", "bce_72h"))) %>%
  mutate(vars = lapply(vars, unlist))
getv <- function(id) MOD$vars[[which(MOD$id == id)]]
PUB <- c(SOFA1 = "SOFA score, 0 h", SOFA3 = "SOFA score, 72 h",
         CI1 = "Complex I-linked respiration, 0 h", BCE1 = "Coupling efficiency, 0 h",
         D1 = "Complex I + coupling, 0 h",
         CI3 = "Complex I-linked respiration, 72 h", BCE3 = "Coupling efficiency, 72 h",
         D3 = "Complex I + coupling, 72 h",
         TCI = "Complex I trajectory", TBCE = "Coupling trajectory",
         TT = "Complex I + coupling trajectories")

# ---- engine -----------------------------------------------------------------
auc_of <- function(y, p) if (length(unique(y)) < 2) NA_real_ else
  as.numeric(pROC::auc(pROC::roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)))
ffit <- function(f, d) tryCatch(logistf::logistf(f, data = d, pl = FALSE,
  control = logistf::logistf.control(maxit = 500, maxstep = 5)), error = function(e) NULL)
pred_f <- function(m, f, d) { X <- model.matrix(f, d); b <- coef(m)[colnames(X)]
  as.numeric(plogis(as.numeric(X %*% b))) }
zv <- function(vs) paste0("z_", vs)

prep <- function(keep, y, extra = character(0)) {
  allv <- unique(c(unlist(MOD$vars), extra))
  d <- df[keep, c(y, allv)]
  for (v in allv) d[[paste0("z_", v)]] <- as.numeric(scale(d[[v]]))
  d
}

cv_oof <- function(vs, d, y) {
  d <- d[complete.cases(d[, c(y, vs)]), ]
  f <- as.formula(paste(y, "~", paste(vs, collapse = "+")))
  n <- nrow(d); A <- numeric(0); P <- matrix(NA_real_, n, R_CV)
  for (r in seq_len(R_CV)) {
    fold <- sample(rep(seq_len(K_CV), length.out = n)); p <- rep(NA_real_, n); bad <- FALSE
    for (k in seq_len(K_CV)) {
      tr <- fold != k
      if (length(unique(d[[y]][tr])) < 2) { bad <- TRUE; break }
      m <- ffit(f, d[tr, , drop = FALSE]); if (is.null(m)) { bad <- TRUE; break }
      pk <- tryCatch(pred_f(m, f, d[!tr, , drop = FALSE]), error = function(e) NULL)
      if (is.null(pk) || any(!is.finite(pk))) { bad <- TRUE; break }
      p[!tr] <- pk
    }
    if (bad || any(is.na(p))) next
    P[, r] <- p; A <- c(A, auc_of(d[[y]], p))
  }
  ok <- colSums(is.na(P)) == 0
  list(area = if (length(A)) mean(A) else NA_real_, n = n, ev = sum(d[[y]]),
       oof = if (any(ok)) rowMeans(P[, ok, drop = FALSE]) else rep(NA_real_, n),
       y = d[[y]])
}

# penalised likelihood ratio test: added coefficients constrained to zero
# within the full Firth model
p_pen <- function(full, added, d, y) {
  d <- d[complete.cases(d[, c(y, full)]), ]
  m <- ffit(as.formula(paste(y, "~", paste(full, collapse = "+"))), d)
  if (is.null(m)) return(NA_real_)
  t <- tryCatch(logistf::logistftest(m, test = as.formula(
         paste("~", paste(added, collapse = "+"), "- 1"))), error = function(e) NULL)
  if (is.null(t)) NA_real_ else as.numeric(t$prob)
}

delong_oof <- function(a, b) {
  if (any(is.na(a$oof)) || any(is.na(b$oof)) || length(a$y) != length(b$y)) return(NA_real_)
  roc_ <- function(z) pROC::roc(z$y, z$oof, levels = c(0, 1), direction = "<", quiet = TRUE)
  t <- tryCatch(suppressWarnings(pROC::roc.test(roc_(a), roc_(b), method = "delong",
                                                paired = TRUE)), error = function(e) NULL)
  if (is.null(t)) NA_real_ else t$p.value
}

rho <- function(x, y) { ok <- complete.cases(x, y)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  c(rho = unname(ct$estimate), P = ct$p.value) }

EV37 <- vapply(ORD, function(y) sum(df[[y]][PAIRED], na.rm = TRUE), numeric(1))
EV40 <- vapply(ORD, function(y) sum(df[[y]], na.rm = TRUE), numeric(1))
HZ37 <- paste0(nice(ORD), " (", EV37, ")")

# =============================================================================
# Table S1A  (no random numbers)
# =============================================================================
L <- df[PAIRED, ]
PNAME <- c(r = "LEAK respiration", ci = "Complex I-linked respiration",
           cii = "Complex II-linked respiration", bce = "Coupling efficiency")
S1A_rtm <- dplyr::bind_rows(lapply(names(PNAME), function(p) {
  a <- L[[paste0(p, "_0h")]]; b <- L[[paste0(p, "_72h")]]
  r1 <- rho(a, b - a); r2 <- rho((a + b) / 2, b - a)
  tibble::tibble(Parameter = unname(PNAME[p]),
         `Baseline-change rho` = f3(r1["rho"]), P = fmt_p(r1["P"]),
         `Oldham rho` = f3(r2["rho"]), `P (Oldham)` = fmt_p(r2["P"])) }))
FLX <- c("r", "ci", "cii")
S1A_flux <- dplyr::bind_rows(lapply(FLX, function(p) {
  r <- rho(L[[paste0("d_", p)]], L$d_bce)
  tibble::tibble(`Change in coupling against` = paste0("change in ", PNAME[p]),
                 rho = f3(r["rho"]), P = fmt_p(r["P"])) }))
fl <- summary(lm(d_bce ~ d_r + d_ci + d_cii, L))
fl_p <- pf(fl$fstatistic[1], fl$fstatistic[2], fl$fstatistic[3], lower.tail = FALSE)
S1A_flux <- dplyr::bind_rows(S1A_flux, tibble::tibble(
  `Change in coupling against` = "changes in all three fluxes together (linear regression)",
  rho = paste0("R2 = ", f3(fl$r.squared)), P = fmt_p(fl_p)))
S1A_anc <- dplyr::bind_rows(lapply(c("ci", "bce"), function(p) {
  row <- tibble::tibble(Parameter = ifelse(p == "ci", "Complex I-linked respiration (pmol/s/mg)",
                                           "Coupling efficiency"))
  for (y in ORD) {
    d <- L[, c(y, paste0(p, "_0h"), paste0(p, "_72h"))]; d <- d[complete.cases(d), ]
    s <- summary(lm(as.formula(paste0(p, "_72h ~ ", p, "_0h + ", y)), d))$coefficients
    row[[nice(y)]] <- paste0(formatC(s[y, 1], format = "f", digits = ifelse(p == "ci", 1, 3)),
                             " (P = ", fmt_p(s[y, 4]), ")")
  }
  row }))

# =============================================================================
# Table S1B  cross-validated AUC, 13 models x 4 horizons
# =============================================================================
FITS <- list()
AR <- dplyr::bind_rows(lapply(ORD, function(y) {
  d <- prep(PAIRED, y)
  dplyr::bind_rows(lapply(seq_len(nrow(MOD)), function(i) {
    s <- cv_oof(zv(MOD$vars[[i]]), d, y)
    FITS[[paste(y, MOD$id[i])]] <<- s
    tibble::tibble(endpoint = nice(y), id = MOD$id[i], area = s$area) })) }))

# admission models on all 40 patients (table S3)
S40a <- dplyr::bind_rows(lapply(ORD, function(y) {
  d <- prep(rep(TRUE, nrow(df)), y)
  dplyr::bind_rows(lapply(c("SOFA1", "CI1", "BCE1", "D1"), function(id) {
    s <- cv_oof(zv(getv(id)), d, y)
    tibble::tibble(endpoint = nice(y), id = id, area40 = s$area) })) }))

# =============================================================================
# Table S1C  against the SOFA score
# =============================================================================
CMP <- list(c("D1","SOFA1"), c("CI1","SOFA1"), c("D3","SOFA1"), c("D3","SOFA3"),
            c("TCI","SOFA1"), c("TCI","SOFA3"), c("TT","SOFA1"), c("TT","SOFA3"),
            c("TBCE","SOFA1"), c("TBCE","SOFA3"), c("D13","SOFA1"))
VS <- dplyr::bind_rows(lapply(ORD, function(y) {
  d <- prep(PAIRED, y)
  dplyr::bind_rows(lapply(CMP, function(cc) {
    a <- FITS[[paste(y, cc[1])]]; b <- FITS[[paste(y, cc[2])]]
    full <- zv(c(getv(cc[2]), getv(cc[1]))); added <- zv(getv(cc[1]))
    sf <- cv_oof(full, d, y)
    tibble::tibble(endpoint = nice(y), m = cc[1], ref = cc[2],
                   area_model = a$area, area_ref = b$area,
                   P_better = delong_oof(a, b), P_adds = p_pen(full, added, d, y),
                   gain = sf$area - b$area) })) }))

# =============================================================================
# Table S1D  within respiration
# =============================================================================
NEST <- list(
  list(full = "D1", add = "bce_0h",  base = "CI1",  q = "Coupling added to Complex I, 0 h"),
  list(full = "D1", add = "ci_0h",   base = "BCE1", q = "Complex I added to coupling, 0 h"),
  list(full = "D3", add = "bce_72h", base = "CI3",  q = "Coupling added to Complex I, 72 h"),
  list(full = "D3", add = "ci_72h",  base = "BCE3", q = "Complex I added to coupling, 72 h"),
  list(full = "TT", add = "rel_bce", base = "TCI",  q = "Coupling trajectory added to Complex I trajectory"),
  list(full = "TT", add = "rel_ci",  base = "TBCE", q = "Complex I trajectory added to coupling trajectory"))
NR <- dplyr::bind_rows(lapply(ORD, function(y) {
  d <- prep(PAIRED, y)
  dplyr::bind_rows(lapply(NEST, function(z) {
    a <- FITS[[paste(y, z$full)]]; b <- FITS[[paste(y, z$base)]]
    tibble::tibble(endpoint = nice(y), q = z$q,
                   P = p_pen(zv(getv(z$full)), zv(z$add), d, y),
                   gain = a$area - b$area) })) }))
NTEST <- nrow(VS) + nrow(NR)
BONF <- 0.05 / NTEST

SHOW1 <- c("SOFA1","SOFA3","CI1","BCE1","D1","CI3","BCE3","D3","TCI","TBCE","TT")
S1B <- AR %>% filter(id %in% SHOW1) %>%
  mutate(Model = unname(PUB[id]), area = f3(area)) %>%
  select(id, Model, endpoint, area) %>%
  pivot_wider(names_from = endpoint, values_from = area) %>%
  dplyr::slice(match(SHOW1, id)) %>% select(-id)
names(S1B)[-1] <- HZ37

SHOWC <- list(c("CI1","SOFA1"), c("D1","SOFA1"), c("D3","SOFA1"), c("D3","SOFA3"),
              c("TBCE","SOFA1"), c("TBCE","SOFA3"), c("TCI","SOFA1"), c("TT","SOFA1"))
S1C <- dplyr::bind_rows(lapply(SHOWC, function(cc) {
  r <- VS %>% filter(m == cc[1], ref == cc[2])
  r <- r[match(nice(ORD), r$endpoint), ]
  tibble::tibble(Model = PUB[cc[1]], Comparator = PUB[cc[2]], Horizon = r$endpoint,
         `AUC, model` = f3(r$area_model), `AUC, comparator` = f3(r$area_ref),
         `Paired P (DeLong)` = fmt_p(r$P_better), `P, added to comparator` = fmt_p(r$P_adds),
         `Gain in AUC` = fg(r$gain)) }))

S1D <- NR %>% mutate(endpoint = factor(endpoint, levels = nice(ORD)),
                     q = factor(q, levels = vapply(NEST, `[[`, "", "q"))) %>%
  arrange(q, endpoint) %>%
  transmute(Comparison = as.character(q), Horizon = as.character(endpoint),
            `Penalized P` = fmt_p(P), `Gain in AUC` = fg(gain))

# =============================================================================
# Table S3  admission models, 37 versus 40 patients (after table S1)
# =============================================================================
S3_EXTRA <- list(
  list(id = "PROF", lab = "Admission respiratory profile (LEAK, CI, CII, coupling)",
       v = c("r_0h", "ci_0h", "cii_0h", "bce_0h")),
  list(id = "LACS", lab = "Arterial lactate + SOFA score, 0 h", v = c("lactate_0h", "sofa_0h")),
  list(id = "MITO", lab = "GDF15 + FGF21, 0 h", v = c("gdf15_0h", "fgf21_0h")))
extra_v <- unique(unlist(lapply(S3_EXTRA, `[[`, "v")))
S3x <- dplyr::bind_rows(lapply(ORD, function(y) {
  d37 <- prep(PAIRED, y, extra_v); d40 <- prep(rep(TRUE, nrow(df)), y, extra_v)
  dplyr::bind_rows(lapply(S3_EXTRA, function(z) {
    s37 <- cv_oof(zv(z$v), d37, y); s40 <- cv_oof(zv(z$v), d40, y)
    tibble::tibble(endpoint = nice(y), id = z$id, lab = z$lab,
                   area37 = s37$area, area40 = s40$area, n40 = s40$n) })) }))
S3a <- AR %>% filter(id %in% c("SOFA1", "CI1", "BCE1", "D1")) %>%
  rename(area37 = area) %>% left_join(S40a, by = c("endpoint", "id")) %>%
  mutate(lab = unname(PUB[id]), n40 = nrow(df))
S3 <- dplyr::bind_rows(S3a, S3x %>% select(endpoint, id, lab, area37, area40, n40))
S3w <- S3 %>%
  mutate(cell = paste0(f3(area37), " / ", f3(area40), " (", fg(area37 - area40), ")")) %>%
  select(id, Model = lab, endpoint, cell) %>%
  pivot_wider(names_from = endpoint, values_from = cell) %>% select(-id)
names(S3w)[-1] <- paste0(nice(ORD), " (", EV37, " / ", EV40, ")")

# =============================================================================
# Table 2  (values of tables S1B and S3, 37 patients)
# =============================================================================
T2_ROWS <- c(CI1  = "Complex I-linked respiration, admission",
             PROF = "Admission respiratory profile (LEAK, Complex I, Complex II, coupling)",
             BCE1 = "Coupling efficiency, admission",
             BCE3 = "Coupling efficiency, 72 h",
             TBCE = "Coupling trajectory",
             SOFA1 = "SOFA score, 0 h",
             SOFA3 = "SOFA score, 72 h")
TAB2 <- dplyr::bind_rows(AR, S3x %>% filter(id == "PROF") %>% select(endpoint, id, area = area37)) %>%
  filter(id %in% names(T2_ROWS)) %>%
  mutate(Marker = unname(T2_ROWS[id]), area = f3(area)) %>%
  select(id, Marker, endpoint, area) %>%
  pivot_wider(names_from = endpoint, values_from = area) %>%
  dplyr::slice(match(names(T2_ROWS), id)) %>% select(Marker, dplyr::all_of(HORIZONS))
names(TAB2)[-1] <- paste0(HORIZONS, " (", EV37[match(HORIZONS, nice(ORD))], ")")
write_out(TAB2, OUT, "Table2.csv")

# =============================================================================
# Table S2  by recovery of coupling efficiency (no random numbers)
# =============================================================================
S2VARS <- list(
  list("Respirometry", "LEAK respiration, 0 h",              "pmol/s/mg", "r_0h"),
  list("Respirometry", "LEAK respiration, 72 h",             "pmol/s/mg", "r_72h"),
  list("Respirometry", "Complex I-linked respiration, 0 h",  "pmol/s/mg", "ci_0h"),
  list("Respirometry", "Complex I-linked respiration, 72 h", "pmol/s/mg", "ci_72h"),
  list("Respirometry", "Complex II-linked respiration, 0 h", "pmol/s/mg", "cii_0h"),
  list("Respirometry", "Complex II-linked respiration, 72 h","pmol/s/mg", "cii_72h"),
  list("Respirometry", "Coupling efficiency, 0 h",           "",          "bce_0h"),
  list("Mitokine",     "GDF15, 0 h",                         "pg/mL",     "gdf15_0h"),
  list("Mitokine",     "GDF15, 72 h",                        "pg/mL",     "gdf15_72h"),
  list("Mitokine",     "FGF21, 0 h",                         "pg/mL",     "fgf21_0h"),
  list("Mitokine",     "FGF21, 72 h",                        "pg/mL",     "fgf21_72h"),
  list("Adipokine",    "Leptin, 0 h",                        "pg/mL",     "leptin_0h"),
  list("Adipokine",    "Ghrelin, 0 h",                       "pg/mL",     "ghrelin_0h"),
  list("Adipokine",    "Leptin to ghrelin ratio, 0 h",       "",          "lepghre_0h"),
  list("Inflammatory", "TNF-\u03b1, 0 h",                    "pg/mL",     "tnfa_0h"),
  list("Inflammatory", "IL-6, 0 h",                          "pg/mL",     "il6_0h"),
  list("Inflammatory", "IL-10, 0 h",                         "pg/mL",     "il10_0h"),
  list("Inflammatory", "IL-10, 72 h",                        "pg/mL",     "il10_72h"),
  list("Inflammatory", "C-reactive protein, 0 h",            "mg/L",      "crp_0h"),
  list("Clinical",     "Arterial lactate, 0 h",              "mmol/L",    "lactate_0h"),
  list("Clinical",     "SOFA score, 0 h",                    "",          "sofa_0h"),
  list("Clinical",     "SAPS 3 score",                       "",          "saps_3"),
  list("Clinical",     "Age",                                "years",     "idade"))
D2 <- df[PAIRED, ]
D2$grp <- ifelse(D2$d_bce > 0, "Recovered", "Not recovered")
nR <- sum(D2$grp == "Recovered"); nL <- sum(D2$grp == "Not recovered")
mq <- function(x, dg) { x <- x[!is.na(x)]; if (!length(x)) return("NA")
  q <- quantile(x, c(.5, .25, .75))
  paste0(formatC(q[1], format = "f", digits = dg), " (", formatC(q[2], format = "f", digits = dg),
         " to ", formatC(q[3], format = "f", digits = dg), ")") }
S2 <- dplyr::bind_rows(lapply(S2VARS, function(z) {
  if (!z[[4]] %in% names(D2)) return(NULL)
  x <- D2[[z[[4]]]]; med <- median(x, na.rm = TRUE)
  dg <- if (is.na(med)) 1 else if (abs(med) < 1) 3 else if (abs(med) < 20) 1 else 0
  tibble::tibble(Domain = z[[1]], Variable = z[[2]], Units = z[[3]], n = sum(!is.na(x)),
         Recovered = mq(x[D2$grp == "Recovered"], dg),
         `Not recovered` = mq(x[D2$grp == "Not recovered"], dg),
         p = suppressWarnings(wilcox.test(x ~ D2$grp, exact = FALSE)$p.value)) }))
S2 <- S2 %>% group_by(Domain) %>% mutate(fdr = p.adjust(p, "BH")) %>% ungroup() %>%
  mutate(P = fmt_p(p), FDR = fmt_p(fdr)) %>% select(-p, -fdr)
names(S2)[names(S2) == "Recovered"] <- paste0("Recovered (n = ", nR, ")")
names(S2)[names(S2) == "Not recovered"] <- paste0("Not recovered (n = ", nL, ")")

# =============================================================================
# Publication format: horizon order, minus signs, rho, powers of ten
# =============================================================================
SUPS <- c("0"="\u2070","1"="\u00b9","2"="\u00b2","3"="\u00b3","4"="\u2074",
          "5"="\u2075","6"="\u2076","7"="\u2077","8"="\u2078","9"="\u2079")
sup <- function(n) paste(SUPS[strsplit(as.character(n), "")[[1]]], collapse = "")
fmt_cell <- function(x) {
  x <- vapply(as.character(x), function(s) {
    for (k in regmatches(s, gregexpr("[0-9.]+e-0*[0-9]+", s))[[1]]) {
      p <- strsplit(k, "e-0*")[[1]]
      s <- sub(k, paste0(p[1], " \u00d7 10\u207b", sup(p[2])), s, fixed = TRUE) }
    s }, character(1), USE.NAMES = FALSE)
  gsub("(^|[ (/=])-([0-9])", "\\1\u2212\\2", x)
}
pub_fmt <- function(d) {
  d <- as.data.frame(d, check.names = FALSE, stringsAsFactors = FALSE)
  hz <- vapply(names(d), function(n) { k <- which(startsWith(n, HORIZONS))
    if (length(k)) k[1] else NA_integer_ }, integer(1))
  if (any(!is.na(hz)))
    d <- d[, c(which(is.na(hz)), which(!is.na(hz))[order(hz[!is.na(hz)])]), drop = FALSE]
  if ("Horizon" %in% names(d)) {
    key <- names(d)[1:(which(names(d) == "Horizon") - 1)]
    blk <- match(do.call(paste, d[key]), unique(do.call(paste, d[key])))
    d <- d[order(blk, match(d$Horizon, HORIZONS)), , drop = FALSE]
  }
  names(d) <- gsub("\\brho\\b", "\u03c1", names(d))
  d[] <- lapply(d, fmt_cell)
  rownames(d) <- NULL
  d
}
S1A_rtm <- pub_fmt(S1A_rtm); S1A_flux <- pub_fmt(S1A_flux); S1A_anc <- pub_fmt(S1A_anc)
S1B <- pub_fmt(S1B); S1C <- pub_fmt(S1C); S1D <- pub_fmt(S1D)
S2 <- pub_fmt(S2); S3w <- pub_fmt(S3w)
S1A_flux[[2]] <- sub("^R2 = ", "R\u00b2 = ", S1A_flux[[2]])

T1 <- "Table S1. Complex I-linked respiration and coupling efficiency, alone and combined, by mortality horizon."
N1A <- paste0("(A) Regression to the mean and divergence by outcome in the ", sum(PAIRED),
  " patients with paired respirometry. Under regression to the mean the correlation between the admission value and its change is negative, and Oldham's correlation (change against the mean of the two measurements) is close to zero. Spearman correlation of the change in coupling with the changes in each respiratory flux, and the variance of the change in coupling explained by the three together. Analysis of covariance: difference in the 72 h value between non-survivors and survivors after adjustment for the admission value.")
N1B <- paste0("(B) Cross-validated areas under the receiver operating characteristic curve (five folds, 20 repetitions, Firth penalized logistic regression fitted inside each training fold, standardized predictors) in the ", sum(PAIRED),
  " patients with paired respirometry. Numbers of deaths in parentheses. Trajectories are relative changes from admission to 72 h.")
N1C <- "(C) Paired P: DeLong test on pooled out-of-fold predicted probabilities (approximate). P, added to comparator: penalized likelihood ratio test for adding the model's variables to the SOFA score, with added coefficients constrained to zero within the full Firth model. Gain: cross-validated area of the SOFA score plus the model minus that of the SOFA score alone."
N1D <- paste0("(D) Penalized likelihood ratio test for adding the stated parameter, and gain in cross-validated area. Across panels C and D, ", NTEST,
  " tests were performed; nominal P values are reported, with a Bonferroni threshold of P < ", signif(BONF, 2), ".")
T2 <- "Table S2. Respirometry, circulating mediators and clinical variables by recovery of coupling efficiency."
N2 <- paste0("Patients with paired respirometry (n = ", sum(PAIRED),
  "). Recovery: any increase in coupling efficiency from admission to 72 h. Median (interquartile range); Mann-Whitney U test; false discovery rate (Benjamini-Hochberg) within each domain.")
T3 <- "Table S3. Admission models evaluated in the patients with paired respirometry and in all patients."
N3 <- paste0("Cross-validated areas under the receiver operating characteristic curve in the ", sum(PAIRED),
  " patients with paired respirometry / in all ", nrow(df),
  " patients (difference). Numbers of deaths in the two sets in the column headers; the three patients who died before 72 h are included only in the second. Same estimator as Table S1.")

# ---- text file --------------------------------------------------------------
LIN <- character(0)
w  <- function(...) LIN <<- c(LIN, paste0(...))
tb <- function(d) LIN <<- c(LIN, capture.output(print(as.data.frame(d), row.names = FALSE, right = FALSE)))
w(T1); w(""); w(N1A); tb(S1A_rtm); w(""); tb(S1A_flux); w(""); tb(S1A_anc)
w(""); w(N1B); tb(S1B); w(""); w(N1C); tb(S1C); w(""); w(N1D); tb(S1D)
w(""); w(T2); w(N2); tb(S2); w(""); w(T3); w(N3); tb(S3w)
writeLines(LIN, file.path(OUT, "Supplementary_Tables.txt"), useBytes = TRUE)

# ---- Word file --------------------------------------------------------------
ft <- function(d) {
  x <- flextable::flextable(as.data.frame(d))
  x <- flextable::font(x, fontname = "Times New Roman", part = "all")
  x <- flextable::fontsize(x, size = 8.5, part = "all")
  x <- flextable::bold(x, part = "header")
  x <- flextable::align(x, j = seq_len(ncol(d))[-1], align = "center", part = "all")
  x <- flextable::padding(x, padding = 1.5, part = "all")
  flextable::autofit(x)
}
par_t <- function(doc, txt, bold = FALSE, size = 10)
  officer::body_add_fpar(doc, officer::fpar(officer::ftext(txt,
    officer::fp_text(font.family = "Times New Roman", font.size = size, bold = bold))))
add_tab <- function(doc, d, note) {
  doc <- flextable::body_add_flextable(doc, ft(d))
  if (!is.null(note)) doc <- par_t(doc, note, size = 8.5)
  par_t(doc, "")
}
doc <- officer::read_docx()
doc <- par_t(doc, T1, bold = TRUE)
doc <- par_t(doc, "A. Regression to the mean and divergence by outcome", bold = TRUE, size = 9.5)
doc <- add_tab(doc, S1A_rtm, NULL); doc <- add_tab(doc, S1A_flux, NULL)
doc <- add_tab(doc, S1A_anc, N1A)
doc <- par_t(doc, "B. Cross-validated discrimination", bold = TRUE, size = 9.5)
doc <- add_tab(doc, S1B, N1B)
doc <- par_t(doc, "C. Comparison with the SOFA score", bold = TRUE, size = 9.5)
doc <- add_tab(doc, S1C, N1C)
doc <- par_t(doc, "D. Within respiration", bold = TRUE, size = 9.5)
doc <- add_tab(doc, S1D, N1D)
doc <- officer::body_add_break(doc)
doc <- par_t(doc, T2, bold = TRUE); doc <- add_tab(doc, S2, N2)
doc <- officer::body_add_break(doc)
doc <- par_t(doc, T3, bold = TRUE); doc <- add_tab(doc, S3w, N3)
doc <- officer::body_end_block_section(doc, officer::block_section(
  officer::prop_section(page_size = officer::page_size(orient = "landscape"), type = "continuous")))
print(doc, target = file.path(OUT, "Supplementary_Tables.docx"))
