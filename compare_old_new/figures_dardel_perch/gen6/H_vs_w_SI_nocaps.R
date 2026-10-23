## SI figure: posterior of harvest intensity H against the cost of reproduction w (final perch ABC).
## Accepted particles sized and shaded by weight, weighted kernel-density HPD contours (50%, 90%), marginal medians.
## Usage: Rscript H_vs_w_SI.R <REEIM_res.rds> <generation> <out_prefix>
suppressMessages({library(data.table); library(ggplot2)})
a <- commandArgs(TRUE); G <- as.integer(a[2]); out <- a[3]
r <- readRDS(a[1]); k <- r$particles[[G]][keep == TRUE][, wn := weights / sum(weights)]
xs <- seq(0.1, 1.5, length.out = 140); ys <- seq(0, 0.2, length.out = 140); ess <- 1 / sum(k$wn^2)
sdw <- function(x) sqrt(sum(k$wn * (x - sum(k$wn * x))^2)); bx <- sdw(k$w) * ess^(-1/6); by <- sdw(k$H) * ess^(-1/6)
z <- outer(xs, ys, Vectorize(function(x, y) sum(k$wn * dnorm(x, k$w, bx) * dnorm(y, k$H, by))))
zs <- sort(as.vector(z), decreasing = TRUE); cz <- cumsum(zs) / sum(zs); lev <- sapply(c(0.5, 0.9), function(p) zs[which(cz >= p)[1]])
dens <- data.table(expand.grid(w = xs, H = ys))[, z := as.vector(z)]
wmed <- function(x) { o <- order(x); x[o][which(cumsum(k$wn[o]) >= 0.5)[1]] }
med <- data.table(w = wmed(k$w), H = wmed(k$H))
p <- ggplot(k[order(wn)], aes(w, H)) +
  geom_point(aes(size = wn, colour = wn), alpha = .8) +
  geom_contour(data = dens, aes(w, H, z = z), breaks = lev, colour = "black", linewidth = 0.8) +
  geom_point(data = med, aes(w, H), shape = 4, size = 5, stroke = 1.6, colour = "#D55E00") +
  scale_colour_viridis_c(name = "posterior\nweight") + scale_size_area(max_size = 4, guide = "none") +
  coord_cartesian(xlim = c(0.1, 1.5), ylim = c(0, 0.2), expand = FALSE) +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank()) +
  labs(x = "cost of reproduction w", y = "harvest intensity H (fraction of biomass removed per year)")
ggsave(paste0(out, ".pdf"), p, width = 7.5, height = 5.5); ggsave(paste0(out, ".png"), p, width = 7.5, height = 5.5, dpi = 200)
