# Tests for discrete subgroups in the four admission respirometry parameters
# (LEAK, Complex I, Complex II, coupling efficiency), standardised, on the raw
# and log1p scales: gap statistic (Ward clustering, 200 reference sets),
# mean silhouette width for two and three clusters, Gaussian mixture model
# (1 to 4 components, BIC), bootstrap stability of the two-cluster solution
# (Jaccard, 500 resamples) and Hartigan's dip test for each parameter.
# Output: output/subgroup_tests/

source("R/00_setup.R")
suppressPackageStartupMessages(library(mclust))   # Mclust needs mclustBIC on the search path
df <- load_data()
OUT <- out_dir("subgroup_tests")
set.seed(42)

RESP <- c("r_0h", "ci_0h", "cii_0h", "bce_0h")
prep <- function(log_scale) {
  d <- na.omit(df[, RESP])
  if (log_scale) d <- as.data.frame(lapply(d, function(x) log1p(pmax(x, 0))))
  scale(as.matrix(d))
}
ward <- function(x, k) list(cluster = cutree(hclust(dist(x), "ward.D2"), k))

structure_tests <- function(X, label) {
  gp <- cluster::clusGap(X, FUN = ward, K.max = 5, B = 200, verbose = FALSE)
  sil <- function(k) mean(cluster::silhouette(ward(X, k)$cluster, dist(X))[, 3])
  mc <- Mclust(X, G = 1:4, verbose = FALSE)
  cb <- fpc::clusterboot(X, B = 500, clustermethod = fpc::hclustCBI, method = "ward.D2",
                         k = 2, seed = 42, count = FALSE)
  tibble::tibble(scale = label, n = nrow(X),
    k_gap = cluster::maxSE(gp$Tab[, "gap"], gp$Tab[, "SE.sim"], method = "Tibs2001SEmax"),
    silhouette_k2 = sil(2), silhouette_k3 = sil(3),
    mixture_components = mc$G, mixture_model = mc$modelName,
    jaccard_cluster1 = cb$bootmean[1], jaccard_cluster2 = cb$bootmean[2])
}
ST <- dplyr::bind_rows(structure_tests(prep(FALSE), "raw"),
                       structure_tests(prep(TRUE),  "log1p"))
write_out(ST, OUT, "structure_tests.csv")

DIP <- dplyr::bind_rows(lapply(RESP, function(v) {
  x <- df[[v]][!is.na(df[[v]])]
  tibble::tibble(parameter = v, n = length(x),
                 dip = unname(diptest::dip.test(x)$statistic),
                 P = diptest::dip.test(x)$p.value,
                 P_log = diptest::dip.test(log1p(pmax(x, 0)))$p.value)
}))
write_out(DIP, OUT, "dip_tests.csv")
