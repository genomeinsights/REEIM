## Monte Carlo error of the final perch posterior summaries (GEN_6). Args: REEIM_res.rds, perch_gen6_traits.txt, output csv.
suppressMessages(library(data.table)); a <- commandArgs(TRUE); set.seed(20261023)
k <- readRDS(a[1])$particles[[6]][keep == TRUE][, wn := weights / sum(weights)]
tr <- fread(a[2], header = FALSE, col.names = c("tf", "LT0", "LT1", "g0", "g1", "h0", "h1"))[, file := sub("^traitStats_", "", tf)]
k <- merge(k, tr, by = "file"); k[, `:=`(dLT = 100 * (LT1 / LT0 - 1), dh = 100 * (h1 / h0 - 1), dg = 100 * (g1 / g0 - 1))]
wq <- function(x, w, p) { o <- order(x); cw <- cumsum(w[o]) / sum(w); x[o][which(cw >= p)[1]] }
st <- list(H_median = c("H", .5), H_q95 = c("H", .95), H_q975 = c("H", .975), w_median = c("w", .5),
           dLT_median = c("dLT", .5), dhmax_median = c("dh", .5), dg_median = c("dg", .5))
est <- sapply(st, function(s) wq(k[[s[1]]], k$wn, as.numeric(s[2])))
B <- 10000; ess <- round(1 / sum(k$wn^2)); n <- nrow(k)
bA <- t(replicate(B, { i <- sample.int(n, ess, replace = TRUE, prob = k$wn); sapply(st, function(s) quantile(k[[s[1]]][i], as.numeric(s[2]), type = 1)) }))
bB <- t(replicate(B, { i <- sample.int(n, n, replace = TRUE); sapply(st, function(s) wq(k[[s[1]]][i], k$wn[i], as.numeric(s[2]))) }))
out <- data.table(statistic = names(st), estimate = est,
  SE_ESS_resample = apply(bA, 2, sd), lo_ESS = apply(bA, 2, quantile, .025), hi_ESS = apply(bA, 2, quantile, .975),
  SE_weighted_bootstrap = apply(bB, 2, sd), lo_wboot = apply(bB, 2, quantile, .025), hi_wboot = apply(bB, 2, quantile, .975))
out[, note := c("", "upper bound of the central 90% interval", "upper bound of the central 95% interval", "", "% change 1900-2025", "% change 1900-2025", "% change 1900-2025")]
out[, `:=`(n_particles = n, ESS = 1 / sum(k$wn^2), B = B)]
fwrite(out, a[3]); print(out[, 1:8], digits = 3)
