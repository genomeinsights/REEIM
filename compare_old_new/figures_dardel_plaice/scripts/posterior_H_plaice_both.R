## Plaice: prior and posterior of harvest intensity H, both final ABC runs, with the ICES assessed range.
suppressMessages({library(data.table); library(ggplot2); library(gridExtra)})
a <- commandArgs(TRUE)   # <res plaice> <res plaice_Hfree> <out prefix>
get <- function(f) { r <- readRDS(f); k <- r$particles[[length(r$particles)]][keep == TRUE]; k[, wn := weights / sum(weights)]; k }
A <- get(a[2]); B <- get(a[1])
th <- theme_bw(base_size = 15) + theme(aspect.ratio = 1, panel.grid.minor = element_blank())
ass <- annotate("rect", xmin = 7.8, xmax = 10.5, ymin = -Inf, ymax = Inf, fill = "#56B4E9", alpha = 0.18)
ass_l <- geom_vline(xintercept = 7.8, colour = "#0072B2", linewidth = 0.6)
pA <- ggplot(A, aes(100 * H)) + ass + ass_l +
  geom_density(aes(weight = wn), fill = "grey", alpha = 0.7, linewidth = 0.6, bounds = c(0, 20)) +
  annotate("segment", x = 0, xend = 20, y = 1 / 20, yend = 1 / 20, linetype = 2, linewidth = 0.6) +
  coord_cartesian(xlim = c(0, 20)) + th +
  labs(x = "Harvest rate (% removed annually)", y = "Density",
       title = "A Prior and posterior of H,\nuninformative prior U(0, 0.2)")
pB <- ggplot(B, aes(100 * H)) + ass + ass_l +
  geom_density(aes(weight = wn), fill = "grey", alpha = 0.7, linewidth = 0.6, bounds = c(6, 11)) +
  annotate("segment", x = 6, xend = 11, y = 1 / 5, yend = 1 / 5, linetype = 2, linewidth = 0.6) +
  coord_cartesian(xlim = c(0, 20)) + th +
  labs(x = "Harvest rate (% removed annually)", y = "Density",
       title = "B Prior and posterior of H,\nconditioned on the assessment U(0.06, 0.11)")
g <- arrangeGrob(pA, pB, nrow = 1)
ggsave(paste0(a[3], ".pdf"), g, width = 12, height = 6); ggsave(paste0(a[3], ".png"), g, width = 12, height = 6, dpi = 150)
