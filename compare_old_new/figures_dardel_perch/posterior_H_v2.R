## Revised posterior_H figure: (A) prior vs posterior of H; (B) posterior vs lake-level empirical
## harvest estimates on the same absolute scale (replaces the z-scored shape comparison).
## Usage (repo root): Rscript compare_old_new/figures_dardel_perch/posterior_H_v2.R [out_folder]
suppressMessages({library(data.table); library(ggplot2); library(gridExtra)})
a <- commandArgs(TRUE); base <- if (length(a) >= 1) a[1] else "./REEIM/out_perch_dardel/"; tag <- if (length(a) >= 2) a[2] else ""; hi <- if (length(a) >= 3) as.numeric(a[3]) else NA
out  <- "./compare_old_new/figures_dardel_perch/"
pr   <- fread("./REEIM/input/perch/priors.csv")[parameter == "H"]
r    <- readRDS(file.path(base, "REEIM_res.rds")); k <- r$particles[[length(r$particles)]][keep == TRUE]
k[, w := weights / sum(weights)]; k[, Hp := 100 * H]
emp  <- fread("./empirical_data/harvestrate.csv")[, .(Hp = harvestperc)]
wq   <- function(x, w, p) { o <- order(x); cw <- cumsum(w[o]); sapply(p, function(q) x[o][which(cw >= q)[1]]) }
ci   <- wq(k$Hp, k$w, c(.05, .5, .95)); inside <- emp[, mean(Hp >= ci[1] & Hp <= ci[3])]
lims <- 100 * c(pr$lower, if (is.na(hi)) pr$upper else hi)
cat(sprintf("posterior median %.1f%%, 90%% CI %.1f-%.1f%%; empirical n=%d, median %.1f%%, %.0f%% inside the 90%% CI, %d above\n",
            ci[2], ci[1], ci[3], nrow(emp), median(emp$Hp), 100 * inside, emp[, sum(Hp > ci[3])]))
th <- theme_bw(base_size = 15) + theme(aspect.ratio = 1, panel.grid.minor = element_blank())
pA <- ggplot(k, aes(Hp)) +
  geom_density(aes(weight = w), fill = "grey", alpha = 0.7, linewidth = 0.6, bounds = lims) +
  annotate("segment", x = lims[1], xend = lims[2], y = 1 / diff(lims), yend = 1 / diff(lims), linetype = 2, linewidth = 0.6) +
  coord_cartesian(xlim = c(0, max(lims[2], max(emp$Hp)))) + th +
  labs(x = "Harvest rate (% removed annually)", y = "Density", title = "A Prior and posterior distributions\nof harvest intensity (H)")
pB <- ggplot() +
  annotate("rect", xmin = ci[1], xmax = ci[3], ymin = -Inf, ymax = Inf, fill = "grey85") +
  geom_density(data = k, aes(Hp, weight = w), linewidth = 0.6, bounds = lims) +
  geom_density(data = emp, aes(Hp), linetype = 2, linewidth = 0.6, bounds = c(0, Inf)) +
  geom_rug(data = emp, aes(Hp), length = unit(0.03, "npc")) +
  geom_vline(xintercept = ci[2], linewidth = 0.4) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5, size = 4.2,
           label = sprintf("%d lakes: median %.1f%%\n%.0f%% within the posterior 90%% CI", nrow(emp), median(emp$Hp), 100 * inside)) +
  coord_cartesian(xlim = c(0, max(lims[2], max(emp$Hp)))) + th +
  labs(x = "Harvest rate (% removed annually)", y = "Density",
       title = "B Empirical harvest rates lie\nwithin the inferred range")
g <- arrangeGrob(pA, pB, nrow = 1)
ggsave(file.path(out, paste0("posterior_H_v2", tag, ".pdf")), g, width = 12, height = 6)
ggsave(file.path(out, paste0("posterior_H_v2", tag, ".png")), g, width = 12, height = 6, dpi = 150)
