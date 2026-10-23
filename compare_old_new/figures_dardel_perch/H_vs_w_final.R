## SI figure: inferred harvest vs cost of reproduction (final perch ABC), and the response per unit harvest by w.
## Usage: Rscript H_vs_w_final.R <REEIM_res.rds> <traits.txt> <generation> <out_prefix>
##   traits.txt: one line per run of that generation: traitStats file name, LT 1900, LT 2025, g 1900, g 2025, hmax 1900, hmax 2025
suppressMessages({library(data.table); library(ggplot2); library(gridExtra)})
a <- commandArgs(TRUE); G <- as.integer(a[3]); out <- a[4]
wq <- function(x, w, p) { o <- order(x); cw <- cumsum(w[o]) / sum(w); sapply(p, function(q) x[o][which(cw >= q)[1]]) }
r <- readRDS(a[1]); k <- r$particles[[G]][keep == TRUE]; k[, wn := weights / sum(weights)]; k[, wg := w * g]
tr <- fread(a[2], header = FALSE, col.names = c("tf", "LT0", "LT1", "g0", "g1", "h0", "h1"))[, file := sub("^traitStats_", "", tf)]
k <- merge(k, tr[, .(file, LT0, LT1)], by = "file"); k[, dLT := 100 * (LT1 / LT0 - 1)]
rHw <- cov.wt(cbind(k$H, k$w), wt = k$wn, cor = TRUE)$cor[1, 2]
br <- seq(0.1, 1.5, length.out = 7)
k[, wbin := cut(w, br, include.lowest = TRUE)]
s <- k[, .(w_mid = wq(w, wn, .5), mass = sum(wn), H_med = wq(H, wn, .5), H_lo = wq(H, wn, .05), H_hi = wq(H, wn, .95)), by = wbin][order(w_mid)]
th <- theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
pA <- ggplot(k, aes(w, H)) +
  geom_vline(xintercept = c(0.6, 1.0), linetype = 2, colour = "grey50") +
  annotate("text", x = c(0.6, 1.0), y = 0.199, label = c("old cap 0.6", "old cap 1.0"), hjust = -0.05, size = 3.2, colour = "grey40") +
  geom_point(aes(size = wn, colour = wg), alpha = .7) + scale_size_area(max_size = 3, guide = "none") +
  scale_colour_viridis_c(name = "extra annual\nmortality of\nmature fish (w x g)") +
  geom_ribbon(data = s, aes(x = w_mid, ymin = H_lo, ymax = H_hi), inherit.aes = FALSE, fill = "grey30", alpha = .15) +
  geom_line(data = s, aes(w_mid, H_med), linewidth = 1) + geom_point(data = s, aes(w_mid, H_med), size = 2) +
  annotate("text", x = 0.12, y = 0.005, hjust = 0, label = sprintf("weighted r(H, w) = %.2f", rHw), size = 3.8) +
  coord_cartesian(ylim = c(0, 0.2)) + th +
  labs(x = "cost of reproduction w", y = "harvest intensity H (fraction of biomass per year)",
       title = sprintf("A  Inferred harvest and cost of reproduction (perch, GEN_%d)", G),
       subtitle = "accepted particles (size = weight); line, band: weighted median and 90% interval per w bin")
tq <- wq(k$w, k$wn, c(1/3, 2/3)); k[, wt3 := factor(fifelse(w < tq[1], "low w", fifelse(w < tq[2], "middle w", "high w")), levels = c("low w", "middle w", "high w"))]
fits <- k[, { m <- lm(dLT ~ H, weights = wn); .(slope = coef(m)[2] / 10, intercept = coef(m)[1]) }, by = wt3]
fits <- fits[order(wt3)]; lab <- paste0(fits$wt3, ": ", sprintf("%.1f", fits$slope), "% per 0.1 H", collapse = "\n")
pB <- ggplot(k, aes(H, dLT, colour = wt3)) + geom_hline(yintercept = 0, linetype = 3) +
  geom_point(aes(size = wn), alpha = .55) + scale_size_area(max_size = 3, guide = "none") +
  geom_smooth(aes(weight = wn), method = "lm", se = FALSE, formula = y ~ x, linewidth = 1) +
  scale_colour_manual(values = c("#0072B2", "#E69F00", "#CC79A7"), name = "cost of reproduction\n(weighted terciles)") +
  annotate("text", x = 0.002, y = min(k$dLT), hjust = 0, vjust = 0, label = lab, size = 3.6) + th +
  labs(x = "harvest intensity H (fraction of biomass per year)", y = "genetic change in length at maturation L_T, 1900-2025 (%)",
       title = "B  Response per unit harvest is similar across w", subtitle = "weighted linear fit per tercile of w")
g <- arrangeGrob(pA, pB, nrow = 1, widths = c(1.15, 1))
ggsave(paste0(out, ".pdf"), g, width = 15, height = 6); ggsave(paste0(out, ".png"), g, width = 15, height = 6, dpi = 150)
cat("weighted r(H, w):", round(rHw, 2), "\nper w bin:\n"); print(s[, .(wbin, mass = round(mass, 3), H_med = round(H_med, 3), H_lo = round(H_lo, 3), H_hi = round(H_hi, 3))])
cat("dLT slope per 0.1 H by w tercile:\n"); print(fits[, .(wt3, slope_pct_per_0.1H = round(slope, 2))])
fwrite(s, paste0(out, "_bins.csv"))
