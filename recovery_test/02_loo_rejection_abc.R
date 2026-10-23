############################################################
##  In-silico recovery test, step 2: leave-one-out rejection ABC
##
##  Each generation-1 prior draw j is treated in turn as pseudo-
##  observed data. Its distance to every other draw i is computed
##  exactly as in run_REEIM_ABC.R,
##     dist = SRMSE_dens + SRMSE_growth + SAD_LT,
##  with draw j's simulated summaries in place of the perch data.
##  The k nearest draws (j itself excluded) form the rejection-ABC
##  posterior for j, as in abc::cv4abc(method = "rejection").
##
##  Reported per quantity and tolerance: bias, prediction error
##  (1 = no better than the prior mean), coverage of the central
##  50/80/90/95% intervals, and Pearson r between true and inferred.
##
##  Run from the repository root:
##    Rscript recovery_test/02_loo_rejection_abc.R
##  Outputs in recovery_test/output/
############################################################

library(data.table)
library(ggplot2)
source("./recovery_test/R/functions.R")

out_dir <- "./recovery_test/output"
gen1    <- readRDS(file.path(out_dir, "gen1_stats.rds"))
load("./empirical_data/observed_data.RData")   # age_structure_obs, size_structure_obs

theta <- gen1$theta
N     <- nrow(theta)

## Tolerances as accepted fractions of the prior draws. 20% is the
## generation-1 cutoff (alpha in run_REEIM_ABC.R); 1% is closer to the
## sharpness of the final SMC posterior.
tol <- c(0.20, 0.05, 0.01)
k   <- pmax(10, round(tol * (N - 1)))

## Quantities to recover. H and the nuisance parameters are inferred
## directly; the r_T are emergent, so their "posterior" is the set of
## realised r_T among the accepted draws.
targets <- c(H = "H", f = "f", w = "w", Shmax = "S[h[max]]", SLT = "S[L[T]]",
             r_hmax = "r[T]~(h)", r_g = "r[T]~(g)", r_LT = "r[T]~(L[T])")
theta[, (grep("^r_", names(theta), value = TRUE)) :=
        lapply(.SD, function(r) 100 * (exp(r) - 1)), .SDcols = patterns("^r_")]  # % change over the century

## ---------------------------------------------------------
## Distance matrix D[j, i]: pseudo-observed j vs simulated i,
## using the driver's distance (R/functions.R). Growth is scored
## over the empirical age composition, with the observed within-age
## variance, and scaled by the sd of the observed sizes-at-age.
## ---------------------------------------------------------
gc     <- growth_constants(age_structure_obs, ncol(gen1$growth))
draw   <- function(j) list(theta = gen1$theta[j], dens = gen1$dens[j, , drop = FALSE],
                           growth = gen1$growth[j, , drop = FALSE])
D <- t(vapply(seq_len(N), function(j) dist_to_target(draw(j), gen1, gc), numeric(N)))
diag(D) <- Inf                       # leave one out

## ---------------------------------------------------------
## Rejection ABC for every pseudo-observed draw
## ---------------------------------------------------------
probs <- c(0.025, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.975)
post <- rbindlist(lapply(seq_along(tol), function(t) {
  rbindlist(lapply(seq_len(N), function(j) {
    acc <- order(D[j, ])[seq_len(k[t])]
    rbindlist(lapply(names(targets), function(p) {
      x <- theta[[p]][acc]
      q <- quantile(x, probs, na.rm = TRUE, names = FALSE)
      data.table(tol = tol[t], j = j, param = p, true = theta[[p]][j],
                 mean = mean(x, na.rm = TRUE), q025 = q[1], q05 = q[2], q10 = q[3],
                 q25 = q[4], median = q[5], q75 = q[6], q90 = q[7], q95 = q[8], q975 = q[9])
    }))
  }))
}))
post <- post[is.finite(true)]
fwrite(post, file.path(out_dir, "loo_posteriors.csv"))

## ---------------------------------------------------------
## Summary metrics
## ---------------------------------------------------------
metrics <- post[, .(
  n           = .N,
  bias        = mean(median - true),
  rel_bias    = mean(median - true) / sd(true),
  pred_error  = sum((median - true)^2) / sum((true - mean(true))^2),
  pearson_r   = cor(true, median),
  cov50       = mean(true >= q25  & true <= q75),
  cov80       = mean(true >= q10  & true <= q90),
  cov90       = mean(true >= q05  & true <= q95),
  cov95       = mean(true >= q025 & true <= q975),
  ci90_width  = mean(q95 - q05) / diff(range(true))
), by = .(tol, param)]
fwrite(metrics, file.path(out_dir, "loo_metrics.csv"))
print(metrics[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])

## Recovery of H within bins of the true H (does accuracy depend on H?)
H_bins <- post[param == "H", .(
  n = .N, bias = mean(median - true), cov90 = mean(true >= q05 & true <= q95),
  ci90_width = mean(q95 - q05)
), by = .(tol, H_bin = cut(true, seq(0, 0.3, 0.05), include.lowest = TRUE))][order(tol, H_bin)]
fwrite(H_bins, file.path(out_dir, "loo_H_by_bin.csv"))

## ---------------------------------------------------------
## Figures
## ---------------------------------------------------------
plot_recovery <- function(dat, pars, file, width, height, ncol) {
  dat <- dat[param %in% pars]
  dat[, lab := factor(targets[param], levels = targets[pars])]
  ann <- metrics[tol == unique(dat$tol) & param %in% pars]
  ann[, lab := factor(targets[param], levels = targets[pars])]
  ann[, txt := sprintf("r = %.2f\ncov90 = %.2f\nbias = %.2g", pearson_r, cov90, bias)]
  set.seed(1)
  bars <- dat[, .SD[sample(.N, min(.N, 300))], by = lab]
  p <- ggplot(dat, aes(true, median)) +
    geom_linerange(data = bars, aes(ymin = q05, ymax = q95), colour = "grey75", linewidth = 0.3) +
    geom_point(size = 0.6, alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "firebrick") +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = txt),
              hjust = -0.05, vjust = 1.1, size = 3, inherit.aes = FALSE) +
    facet_wrap(~ lab, scales = "free", labeller = label_parsed, ncol = ncol) +
    labs(x = "True value (pseudo-observed draw)", y = "Inferred (posterior median, 90% interval)") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), strip.background = element_blank())
  ggsave(file, p, width = width, height = height)
}

main_tol <- 0.05
plot_recovery(post[tol == main_tol], c("H", "r_hmax", "r_g", "r_LT"),
              file.path(out_dir, "Fig_S_recovery_H_rT.pdf"), 12, 3.6, 4)
plot_recovery(post[tol == main_tol], c("f", "w", "Shmax", "SLT"),
              file.path(out_dir, "Fig_S_recovery_nuisance.pdf"), 12, 3.6, 4)

## Coverage calibration across tolerances
cov_long <- melt(metrics, id.vars = c("tol", "param"),
                 measure.vars = c("cov50", "cov80", "cov90", "cov95"),
                 variable.name = "nominal", value.name = "empirical")
cov_long[, nominal := as.numeric(sub("cov", "", nominal)) / 100]
cov_long[, lab := factor(targets[param], levels = targets)]
p_cov <- ggplot(cov_long, aes(nominal, empirical, colour = factor(tol))) +
  geom_abline(linetype = 2, colour = "grey50") +
  geom_point() + geom_line() +
  facet_wrap(~ lab, labeller = label_parsed, ncol = 4) +
  coord_equal(xlim = c(0.45, 1), ylim = c(0.45, 1)) +
  labs(x = "Nominal coverage", y = "Empirical coverage", colour = "Accepted\nfraction") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank())
ggsave(file.path(out_dir, "Fig_S_recovery_coverage.pdf"), p_cov, width = 10, height = 5.5)

## LaTeX table for the SI (tolerance used in the main figure)
tab <- metrics[tol == main_tol]
tab[, Quantity := targets[param]]
tex <- c("\\begin{tabular}{lrrrrr}", "\\hline",
         "Quantity & Bias & Rel. bias & Pred. error & Pearson $r$ & 90\\% coverage \\\\", "\\hline",
         tab[, sprintf("%s & %.3g & %.2f & %.2f & %.2f & %.2f \\\\",
                       param, bias, rel_bias, pred_error, pearson_r, cov90)],
         "\\hline", "\\end{tabular}")
writeLines(tex, file.path(out_dir, "Table_S_recovery.tex"))
