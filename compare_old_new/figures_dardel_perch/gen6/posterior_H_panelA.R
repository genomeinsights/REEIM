## Perch: prior and posterior of harvest intensity H (final generation, weighted). Single panel.
suppressMessages({library(data.table); library(ggplot2)})
a <- commandArgs(TRUE); r <- readRDS(a[1]); out <- a[2]
pr <- fread("./REEIM/input/perch/priors.csv")[parameter == "H"]; lims <- 100 * c(pr$lower, pr$upper)
k <- r$particles[[length(r$particles)]][keep == TRUE][, `:=`(wn = weights / sum(weights), Hp = 100 * H)]
p <- ggplot(k, aes(Hp)) +
  geom_density(aes(weight = wn), fill = "grey", alpha = 0.7, linewidth = 0.6, bounds = lims) +
  annotate("segment", x = lims[1], xend = lims[2], y = 1 / diff(lims), yend = 1 / diff(lims), linetype = 2, linewidth = 0.6) +
  coord_cartesian(xlim = lims) + theme_bw(base_size = 15) + theme(aspect.ratio = 1, panel.grid.minor = element_blank()) +
  labs(x = "Harvest rate (% of biomass removed annually)", y = "Density")
ggsave(paste0(out, ".pdf"), p, width = 6, height = 6); ggsave(paste0(out, ".png"), p, width = 6, height = 6, dpi = 200)
