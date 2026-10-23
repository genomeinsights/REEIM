############################################################
##  Summarise the ABC-SMC runs on pseudo-observed data
##  (03_smc_pseudo_obs.R): weighted posterior median and 90%
##  interval per generation, against the known true values.
##
##  Run from the repository root:
##    Rscript recovery_test/04_smc_summary.R
############################################################

library(data.table)
library(ggplot2)

smc_root <- "./recovery_test/smc"
out_dir  <- "./recovery_test/output"
targets  <- c(H = "H", r_hmax = "r[T]~(h)", r_g = "r[T]~(g)", r_LT = "r[T]~(L[T])")

wq <- function(x, w, p) { o <- order(x); x[o][which(cumsum(w[o]) / sum(w) >= p)[1]] }
pct <- function(r) 100 * (exp(r) - 1)

runs <- list.files(smc_root, pattern = "^smc_res.rds$", recursive = TRUE, full.names = TRUE)
summ <- rbindlist(lapply(runs, function(f) {
  smc <- readRDS(f)
  id  <- sub("PO_", "", basename(dirname(f)))
  rbindlist(lapply(seq_along(smc$particles), function(G) {
    p <- smc$particles[[G]][keep == TRUE]
    rbindlist(lapply(names(targets), function(v) {
      x <- p[[v]]; tr <- smc$truth[[v]]
      if (v != "H") { x <- pct(x); tr <- pct(tr) }
      data.table(id = id, gen = G, param = v, true = tr,
                 median = wq(x, p$weights, 0.5),
                 q05 = wq(x, p$weights, 0.05), q95 = wq(x, p$weights, 0.95),
                 ESS = smc$ESS[G])
    }))
  }))
}))
summ[, covered := true >= q05 & true <= q95]
fwrite(summ, file.path(out_dir, "smc_pseudo_obs_summary.csv"))
print(summ[gen == max(gen)])

last <- summ[, .SD[gen == max(gen)], by = id]
last[, lab := factor(targets[param], levels = targets)]
p <- ggplot(last, aes(true, median)) +
  geom_abline(linetype = 2, colour = "firebrick") +
  geom_errorbar(aes(ymin = q05, ymax = q95), width = 0) +
  geom_point() +
  facet_wrap(~ lab, scales = "free", labeller = label_parsed, nrow = 1) +
  labs(x = "True value (pseudo-observed dataset)",
       y = "ABC-SMC posterior (median, 90% interval)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank())
ggsave(file.path(out_dir, "Fig_S_smc_pseudo_obs.pdf"), p, width = 12, height = 3.6)
