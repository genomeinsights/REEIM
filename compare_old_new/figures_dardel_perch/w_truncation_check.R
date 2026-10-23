## Does the w (0.6) / SLT (250) prior truncation pull the upper bound on H (and on dL_T) down?
suppressMessages(library(data.table))
wq <- function(x, w, p) { o <- order(x); cw <- cumsum(w[o]) / sum(w); sapply(p, function(q) x[o][which(cw >= q)[1]]) }
wcor <- function(x, y, w) { mx <- sum(w*x); my <- sum(w*y); sum(w*(x-mx)*(y-my)) / sqrt(sum(w*(x-mx)^2) * sum(w*(y-my)^2)) }
chk <- function(lab, dir) {
  r <- readRDS(file.path(dir, "REEIM_res.rds")); G <- length(r$particles); k <- r$particles[[G]][keep == TRUE]; k[, wn := weights / sum(weights)]
  k[, dLT := sapply(file, function(f) { t <- fread(file.path(dir, paste0("GEN_", G), paste0("traitStats_", f)))[G %in% c(10000, 10125)]; 100 * (t$LT[2] / t$LT[1] - 1) })]
  cat("\n==", lab, " (GEN", G, ", ESS", round(1 / sum(k$wn^2)), ")\n")
  cat("weighted cor(H, .):", paste(sapply(c("w", "SLT", "f", "Shmax"), function(v) sprintf("%s %.2f", v, wcor(k$H, k[[v]], k$wn))), collapse = "  "), "\n")
  bands <- list(`w<0.4` = k$w < 0.4, `0.4-0.5` = k$w >= 0.4 & k$w < 0.5, `0.5-0.6` = k$w >= 0.5, `0.57-0.6` = k$w >= 0.57)
  for (b in names(bands)) { s <- k[bands[[b]]]
    cat(sprintf("  w %-8s mass %.2f  H q50 %.3f q90 %.3f q95 %.3f | dLT q05 %.1f%%\n", b, sum(s$wn), wq(s$H, s$wn, .5), wq(s$H, s$wn, .9), wq(s$H, s$wn, .95), wq(s$dLT, s$wn, .05))) }
  hi <- k$H > wq(k$H, k$wn, .9)
  cat(sprintf("  top-10%% H particles: weighted mean w %.3f vs others %.3f; share with w>=0.57: %.2f vs %.2f\n",
      sum(k$wn[hi] * k$w[hi]) / sum(k$wn[hi]), sum(k$wn[!hi] * k$w[!hi]) / sum(k$wn[!hi]), sum(k$wn[hi & k$w >= .57]) / sum(k$wn[hi]), sum(k$wn[!hi & k$w >= .57]) / sum(k$wn[!hi])))
  cat(sprintf("  SLT>=237.5 (top 5%%): mass %.2f, H q95 %.3f (all: %.3f)\n", k[SLT >= 237.5, sum(wn)], k[SLT >= 237.5, wq(H, wn, .95)], wq(k$H, k$wn, .95)))
  ## weighted linear regression of H on w, as a rough extrapolation of the 95% H bound's dependence on w
  m <- lm(H ~ w, data = k, weights = wn); q <- sapply(seq(0.1, 0.6, 0.1), function(a) { s <- k[abs(w - a) < 0.075]; if (sum(s$wn) > 0.02) wq(s$H, s$wn, .95) else NA })
  cat("  H q95 in w windows (±0.075) at w = 0.1..0.6:", paste(round(q, 3), collapse = " "), "\n")
  cat(sprintf("  weighted lm H ~ w: slope %.3f per unit w (SE %.3f)\n", coef(m)[2], summary(m)$coef[2, 2]))
}
chk("NEW", "REEIM/out_perch_dardel"); chk("OLD", "REEIM/out_REEIM")
