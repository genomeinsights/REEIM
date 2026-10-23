## w-conditional table for the final perch ABC (one generation): by w band, posterior mass, H median and 95%,
## change in L_T 1900-2025 median and 5% quantile, and realised reproduction-linked mortality w*g.
## Usage: Rscript w_conditional.R <REEIM_res.rds> <traits.txt> <generation> <out.csv>
suppressMessages(library(data.table))
a <- commandArgs(TRUE); G <- as.integer(a[3])
wq <- function(x, w, p) { o <- order(x); cw <- cumsum(w[o]) / sum(w); sapply(p, function(q) x[o][which(cw >= q)[1]]) }
r <- readRDS(a[1]); k <- r$particles[[G]][keep == TRUE][, wn := weights / sum(weights)]
tr <- fread(a[2], header = FALSE, col.names = c("tf", "LT0", "LT1", "g0", "g1", "h0", "h1"))[, file := sub("^traitStats_", "", tf)]
k <- merge(k, tr[, .(file, LT0, LT1)], by = "file")[, `:=`(dLT = 100 * (LT1 / LT0 - 1), wg = w * g)]
br <- c(0.1, 0.5, 0.9, 1.2, 1.5)
k[, band := cut(w, br, include.lowest = TRUE, labels = c("0.1-0.5", "0.5-0.9", "0.9-1.2", "1.2-1.5"))]
one <- function(s, lab) data.table(w_band = lab, mass = sum(s$wn), H_median = wq(s$H, s$wn, .5), H_q95 = wq(s$H, s$wn, .95),
  dLT_median = wq(s$dLT, s$wn, .5), dLT_q05 = wq(s$dLT, s$wn, .05), wg_median = wq(s$wg, s$wn, .5))
tab <- rbind(rbindlist(lapply(levels(k$band), function(b) one(k[band == b], b))), one(k, "all"))
print(tab, digits = 3); fwrite(tab, a[4])
