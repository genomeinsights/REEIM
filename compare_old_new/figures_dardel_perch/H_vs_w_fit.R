## H vs w with fit information (final perch ABC).
## A: every valid run of GEN_1..GEN_g coloured by total ABC distance (best on top); the accepted particles of GEN_g ringed.
## B: accepted particles of GEN_g sized and shaded by posterior weight, with weighted 2-D kernel density HPD contours (50%, 90%).
## Usage: Rscript H_vs_w_fit.R <REEIM_res.rds> <generation> <out_prefix>
suppressMessages({library(data.table); library(ggplot2); library(gridExtra)})
a <- commandArgs(TRUE); G <- as.integer(a[2]); out <- a[3]
r <- readRDS(a[1]); eps <- r$epsilons[[G]]
allp <- rbindlist(lapply(seq_len(G), function(g) r$particles[[g]][, .(w, H, dist, gen = g)]))
k <- r$particles[[G]][keep == TRUE][, wn := weights / sum(weights)]
## weighted Gaussian KDE on a grid (bandwidths: weighted Scott's rule), HPD levels from the grid
xs <- seq(0.1, 1.5, length.out = 120); ys <- seq(0, 0.2, length.out = 120); ess <- 1 / sum(k$wn^2)
sdw <- function(x) sqrt(sum(k$wn * (x - sum(k$wn * x))^2)); bx <- sdw(k$w) * ess^(-1/6); by <- sdw(k$H) * ess^(-1/6)
z <- outer(xs, ys, Vectorize(function(x, y) sum(k$wn * dnorm(x, k$w, bx) * dnorm(y, k$H, by))))
zs <- sort(as.vector(z), decreasing = TRUE); cz <- cumsum(zs) / sum(zs)
lev <- sapply(c(0.5, 0.9), function(p) zs[which(cz >= p)[1]])
dens <- data.table(expand.grid(w = xs, H = ys))[, z := as.vector(z)]
wmed <- function(x) { o <- order(x); x[o][which(cumsum(k$wn[o]) >= 0.5)[1]] }
med <- data.table(w = wmed(k$w), H = wmed(k$H))
th <- theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
common <- list(geom_vline(xintercept = c(0.6, 1.0), linetype = 2, colour = "grey50"),
               coord_cartesian(xlim = c(0.1, 1.5), ylim = c(0, 0.2)),
               labs(x = "cost of reproduction w", y = "harvest intensity H (fraction of biomass per year)"), th)
pA <- ggplot(allp[order(-dist)], aes(w, H, colour = dist)) + geom_point(size = 0.7, alpha = .8) +
  scale_colour_viridis_c(direction = -1, name = "ABC distance\n(lower = better)", limits = quantile(allp$dist, c(0, .9)), oob = scales::squish) +
  geom_point(data = k, inherit.aes = FALSE, aes(w, H), colour = "black", shape = 1, size = 1.3, stroke = .3) + common +
  ggtitle(sprintf("A  All %d valid runs of GEN_1-%d, coloured by distance", nrow(allp), G),
          subtitle = sprintf("rings: %d accepted in GEN_%d (distance < %.3f); dashed: old w caps", nrow(k), G, eps))
pB <- ggplot(k[order(wn)], aes(w, H)) + geom_point(aes(size = wn, colour = wn), alpha = .8) +
  geom_contour(data = dens, aes(w, H, z = z), breaks = lev, colour = "black", linewidth = c(0.9)) +
  geom_point(data = med, aes(w, H), inherit.aes = FALSE, shape = 4, size = 5, stroke = 1.6, colour = "red") +
  scale_colour_viridis_c(name = "posterior\nweight") + scale_size_area(max_size = 4, guide = "none") + common +
  ggtitle(sprintf("B  GEN_%d posterior: particles by weight (ESS %.0f)", G, ess),
          subtitle = "contours: weighted density, 50% and 90% HPD; red cross: marginal medians")
g <- arrangeGrob(pA, pB, nrow = 1)
ggsave(paste0(out, ".pdf"), g, width = 15, height = 6); ggsave(paste0(out, ".png"), g, width = 15, height = 6, dpi = 150)
best <- allp[order(dist)][1:500]
cat(sprintf("best 500 of all runs: w median %.2f, H median %.3f, spearman r(H, w) %.2f | all valid runs: r(dist, H) %.2f, r(dist, w) %.2f\n",
    median(best$w), median(best$H), cor(best$H, best$w, method = "spearman"), cor(allp$dist, allp$H, method = "spearman"), cor(allp$dist, allp$w, method = "spearman")))
mx <- dens[which.max(z)]; cat(sprintf("posterior density peak at w %.2f, H %.3f\n", mx$w, mx$H))
