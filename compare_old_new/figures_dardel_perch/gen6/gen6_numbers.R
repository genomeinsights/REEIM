suppressMessages({library(data.table); library(loo)})
a <- commandArgs(TRUE); R6 <- a[1]; tr6 <- a[2]; R5 <- a[3]; tr5 <- a[4]; A <- a[5]
wq <- function(x, w, p) { o <- order(x); cw <- cumsum(w[o]) / sum(w); sapply(p, function(q) x[o][which(cw >= q)[1]]) }
wcor <- function(x, y, w) cov.wt(cbind(x, y), wt = w, cor = TRUE)$cor[1, 2]
r6 <- readRDS(R6); r5 <- readRDS(R5)
## tolerance and ESS per generation
tol <- rbindlist(lapply(seq_along(r6$particles), function(g) { k <- r6$particles[[g]][keep == TRUE]; w <- k$weights / sum(k$weights)
  data.table(generation = g, epsilon = r6$epsilons[[g]], n_accepted = nrow(k), ESS = 1 / sum(w^2), ESS_frac = 1 / sum(w^2) / nrow(k), max_weight = max(w)) }))
fwrite(tol, file.path(A, "tolerance_ESS_per_generation.csv")); print(tol, digits = 3)
load_k <- function(r, trf) { G <- length(r$particles); k <- r$particles[[G]][keep == TRUE][, wn := weights / sum(weights)]
  tr <- fread(trf, header = FALSE, col.names = c("tf", "LT0", "LT1", "g0", "g1", "h0", "h1"))[, file := sub("^traitStats_", "", tf)]
  k <- merge(k, tr, by = "file"); k[, `:=`(dLT = 100 * (LT1 / LT0 - 1), dh = 100 * (h1 / h0 - 1), dg_abs = g1 - g0, dg = 100 * (g1 / g0 - 1), wg = w * g)]; k }
k <- load_k(r6, tr6); k5 <- load_k(r5, tr5)
cat(sprintf("\nGEN_6: ESS %.0f | weighted r(H, w) = %.2f | r(H, w*g) = %.2f | P(H<0.01) = %.3f | P(dLT<0) = %.3f\n", 1 / sum(k$wn^2), wcor(k$H, k$w, k$wn), wcor(k$H, k$wg, k$wn), sum(k$wn[k$H < 0.01]), sum(k$wn[k$dLT < 0])))
for (v in c("H", "w", "f", "SLT", "Shmax", "g", "wg", "dLT", "dh", "dg", "dg_abs")) cat(sprintf("  %-6s %.4g [%.4g, %.4g]\n", v, wq(k[[v]], k$wn, .5), wq(k[[v]], k$wn, .05), wq(k[[v]], k$wn, .95)))
cat(sprintf("  H 97.5%% quantile %.3f; mass H > 0.19: %.3f; mass w > 1.43: %.3f\n", wq(k$H, k$wn, .975), sum(k$wn[k$H > 0.19]), sum(k$wn[k$w > 1.43])))
## per 0.1 H (weighted linear regressions on H)
per <- rbindlist(lapply(c("dLT", "dh", "dg_abs", "dg"), function(v) { m <- lm(k[[v]] ~ k$H, weights = k$wn)
  data.table(quantity = v, per_0.1_H = coef(m)[2] / 10, se = summary(m)$coef[2, 2] / 10, weighted_r = wcor(k$H, k[[v]], k$wn)) }))
fwrite(per, file.path(A, "change_per_0.1H.csv")); cat("\nper 0.1 H:\n"); print(per, digits = 3)
## weight robustness on GEN_6
w0 <- k$wn; N <- length(w0); wt <- pmin(w0, mean(w0) * sqrt(N)); wt <- wt / sum(wt)
ps <- suppressWarnings(psis(log(w0), r_eff = NA)); wp <- as.vector(weights(ps, log = FALSE, normalize = TRUE)); kh <- pareto_k_values(ps)
row <- function(lab, w) data.table(method = lab, ESS = 1 / sum(w^2), max_weight = max(w), H_median = wq(k$H, w, .5), H_q05 = wq(k$H, w, .05), H_q95 = wq(k$H, w, .95),
  mass_H_gt_0.19 = sum(w[k$H > 0.19]), w_median = wq(k$w, w, .5), w_q05 = wq(k$w, w, .05), w_q95 = wq(k$w, w, .95), mass_w_gt_1.43 = sum(w[k$w > 1.43]),
  dLT_median = wq(k$dLT, w, .5), dLT_q05 = wq(k$dLT, w, .05), dLT_q95 = wq(k$dLT, w, .95))
rob <- rbind(row("raw importance weights", w0), row("truncated (Ionides 2008)", wt), row(sprintf("PSIS (k-hat %.2f)", kh), wp))
fwrite(rob, file.path(A, "gen6_weight_robustness.csv")); cat("\nweight robustness:\n"); print(rob, digits = 3)
## lakes inside the posterior 90% interval of H
emp <- fread("./empirical_data/harvestrate.csv")$harvestperc / 100; ci <- wq(k$H, k$wn, c(.05, .95))
cat(sprintf("\nlakes: %d; inside posterior 90%% interval [%.3f, %.3f]: %.0f%% (%d); below: %d; above: %d\n", length(emp), ci[1], ci[2], 100 * mean(emp >= ci[1] & emp <= ci[2]), sum(emp >= ci[1] & emp <= ci[2]), sum(emp < ci[1]), sum(emp > ci[2])))
## GEN_5 vs GEN_6 conditional on H (kernel h = 0.02, as RUN_PPS)
cmp <- rbindlist(lapply(list(c("GEN_5", "k5"), c("GEN_6", "k")), function(z) { kk <- get(z[2])
  rbindlist(lapply(c(0.02, 0.10, 0.16), function(Ht) { ww <- kk$wn * exp(-0.5 * ((kk$H - Ht) / 0.02)^2); ww <- ww / sum(ww)
    data.table(gen = z[1], H = Ht, ESS_cond = 1 / sum(ww^2), w = wq(kk$w, ww, .5), f = wq(kk$f, ww, .5), SLT = wq(kk$SLT, ww, .5), Shmax = wq(kk$Shmax, ww, .5), dLT = wq(kk$dLT, ww, .5), dLT_q05 = wq(kk$dLT, ww, .05), dLT_q95 = wq(kk$dLT, ww, .95)) })) }))
fwrite(cmp, file.path(A, "gen5_vs_gen6_conditional_on_H.csv")); cat("\nGEN_5 vs GEN_6, conditional on H:\n"); print(cmp[order(H, gen)], digits = 3)
