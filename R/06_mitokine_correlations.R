# Cross-sectional Spearman correlations of the circulating mediators GDF15,
# FGF21, leptin and ghrelin with the four respirometry parameters, at 0 h and
# at 72 h (4 x 4 x 2 = 32 tests), false discovery rate (Benjamini-Hochberg)
# across all 32.
# Output: output/mitokine_correlations/correlations.csv

source("R/00_setup.R")
df <- load_data()
OUT <- out_dir("mitokine_correlations")

MED  <- c("gdf15", "fgf21", "leptin", "ghrelin")
RESP <- c("r", "ci", "cii", "bce")
COR <- dplyr::bind_rows(lapply(c("0h", "72h"), function(t)
  dplyr::bind_rows(lapply(MED, function(m) dplyr::bind_rows(lapply(RESP, function(r) {
    x <- df[[paste0(m, "_", t)]]; y <- df[[paste0(r, "_", t)]]; ok <- complete.cases(x, y)
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
    tibble::tibble(time = t, mediator = m, respirometry = r, n = sum(ok),
                   rho = unname(ct$estimate), P = ct$p.value) }))))))
COR <- COR %>% mutate(FDR = p.adjust(P, "BH")) %>% arrange(P)
write_out(COR, OUT, "correlations.csv")
message(sum(COR$FDR < .05), " of ", nrow(COR), " correlations with FDR < 0.05")
