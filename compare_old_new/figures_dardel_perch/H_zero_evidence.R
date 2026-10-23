## Evidence for H = 0 in the perch ABC posterior (weighted final-generation particles).
suppressMessages(library(data.table)); set.seed(1)
wq <- function(x, w, p) { o <- order(x); cw <- cumsum(w[o]) / sum(w); sapply(p, function(q) x[o][which(cw >= q)[1]]) }
## reflection KDE density at 0: f(0) = 2 * sum(w * dnorm(x/bw)) / bw   (mirror at the boundary)
f0 <- function(x, w, bw) 2 * sum(w * dnorm(x / bw)) / bw
bw_nrd <- function(x, w, ess) { m <- sum(w * x); s <- sqrt(sum(w * (x - m)^2)); 0.9 * min(s, (wq(x, w, .75) - wq(x, w, .25)) / 1.34) * ess^-0.2 }
run <- function(lab, file, hi) {
  r <- readRDS(file); k <- r$particles[[length(r$particles)]][keep == TRUE]; x <- k$H; w <- k$weights / sum(k$weights); ess <- 1 / sum(w^2)
  cat(sprintf("\n%s: GEN %d, kept %d, ESS %.0f, H prior U(0, %g)\n", lab, length(r$particles), nrow(k), ess, hi))
  for (d in c(0.005, 0.01, 0.02, 0.03)) cat(sprintf("  P(H<%.3f): posterior %.3f  prior %.3f\n", d, sum(w[x < d]), d / hi))
  q <- wq(x, w, c(.01, .025, .05, .10)); cat(sprintf("  weighted quantiles 1%% %.4f  2.5%% %.4f  5%% %.4f  10%% %.4f\n", q[1], q[2], q[3], q[4]))
  bw0 <- bw_nrd(x, w, ess)
  for (m in c(0.5, 1, 2)) { bw <- m * bw0
    bf <- f0(x, w, bw) * hi
    bs <- replicate(1000, { i <- sample(length(x), round(ess), TRUE, prob = w); f0(x[i], rep(1 / length(i), length(i)), bw) * hi })
    cat(sprintf("  Savage-Dickey BF01 (H=0 vs H>0), reflection KDE bw %.4f (x%.1f): %.2f  [bootstrap 90%%: %.2f-%.2f]\n", bw, m, bf, quantile(bs, .05), quantile(bs, .95))) }
  po <- sum(w[x < 0.01]); pr <- 0.01 / hi
  bs <- replicate(1000, { i <- sample(length(x), round(ess), TRUE, prob = w); p <- mean(x[i] < 0.01); (p / (1 - p)) / (pr / (1 - pr)) })
  cat(sprintf("  BF(H<0.01 vs H>=0.01) = %.2f  [bootstrap 90%%: %.2f-%.2f]\n", (po / (1 - po)) / (pr / (1 - pr)), quantile(bs, .05), quantile(bs, .95)))
}
run("NEW (Dardel)", "REEIM/out_perch_dardel/REEIM_res.rds", 0.2)
run("OLD (first submission)", "REEIM/out_REEIM/REEIM_res.rds", 0.3)
