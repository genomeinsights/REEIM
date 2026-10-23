## Change in L_T (1900-2025) against w, coloured by H; partial effects of H and w.
suppressMessages({library(data.table); library(ggplot2)})
a <- commandArgs(TRUE); r <- readRDS(a[1]); out <- a[3]
k <- r$particles[[length(r$particles)]][keep == TRUE][, wn := weights / sum(weights)]
tr <- fread(a[2], header = FALSE, col.names = c("tf", "LT0", "LT1", "g0", "g1", "h0", "h1"))[, file := sub("^traitStats_", "", tf)]
k <- merge(k, tr, by = "file")[, `:=`(dLT = 100 * (LT1 / LT0 - 1), dh = 100 * (h1 / h0 - 1), dg = 100 * (g1 / g0 - 1))]
wc <- function(x, y) cov.wt(cbind(x, y), wt = k$wn, cor = TRUE)$cor[1, 2]
cat(sprintf("weighted r: dLT~w %.2f, dLT~H %.2f, dh~w %.2f, dg~w %.2f\n", wc(k$dLT, k$w), wc(k$dLT, k$H), wc(k$dh, k$w), wc(k$dg, k$w)))
z <- function(x) (x - sum(k$wn * x)) / sqrt(sum(k$wn * (x - sum(k$wn * x))^2))
for (v in c("dLT", "dh", "dg")) { m <- lm(z(k[[v]]) ~ z(k$H) + z(k$w), weights = k$wn); cat(sprintf("%s ~ H + w (standardised): H %.2f, w %.2f, R2 %.2f\n", v, coef(m)[2], coef(m)[3], summary(m)$r.squared)) }
p <- ggplot(k[order(H)], aes(w, dLT, colour = H)) + geom_hline(yintercept = 0, linetype = 3) + geom_point(aes(size = wn), alpha = .8) +
  scale_size_area(max_size = 3, guide = "none") + scale_colour_viridis_c(name = "harvest\nintensity H") +
  geom_smooth(aes(weight = wn), method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = .8) +
  theme_bw(base_size = 13) + theme(panel.grid.minor = element_blank()) +
  labs(x = "cost of reproduction w", y = "genetic change in length at maturation L_T, 1900-2025 (%)")
ggsave(paste0(out, ".png"), p, width = 7, height = 5, dpi = 180); ggsave(paste0(out, ".pdf"), p, width = 7, height = 5)
