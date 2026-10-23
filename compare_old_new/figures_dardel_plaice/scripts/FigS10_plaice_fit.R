## Fig. S10: fit of REEIM to present-day plaice data (final assessment-conditioned posterior, GEN_5 of out_plaice_final).
## Style as the perch data-fit figure (R/PPC_ABC_PPS.R): accepted particles' simulated size distribution and size at age
## at the comparison cycle (2020) against the observed data, plus maturation (posterior L_T vs observed female L50).
## Envelopes are weighted by the ABC importance weights. Args: REEIM_res.rds, folder with the accepted demographicData files, output prefix.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
a <- commandArgs(TRUE); r <- readRDS(a[1]); dd <- a[2]; out <- a[3]
load("REEIM/plaice_validation/data/observed_data_plaice_female.RData")   # age_structure_obs, size_structure_obs, ogive
APOT <- 10000 + nrow(fread("REEIM/input/plaice/fishing_intensity.csv")) - 1  # 2020
BW <- 2; from <- 16; to <- 50   # input/plaice/abc_settings.csv
k <- r$particles[[5]][keep == TRUE][, wn := weights / sum(weights)]
wq <- function(x, w, p) { o <- order(x); x <- x[o]; cw <- cumsum(w[o]) / sum(w); sapply(p, function(q) x[which(cw >= q)[1]]) }
sims <- rbindlist(lapply(seq_len(nrow(k)), function(i) {
  s <- fread(file.path(dd, paste0("demographicData_", k$file[i])), fill = Inf)
  cs <- as.numeric(s[V1 == "countSize" & V2 == APOT, 3:ncol(s)][1, ]); cs[is.na(cs)] <- 0
  sz <- rep(seq_along(cs) - 1, cs); sz <- sz[sz >= from & sz <= to]; de <- density(sz, bw = BW, from = from, to = to)
  ms <- as.numeric(s[V1 == "meanSize" & V2 == APOT, 3:ncol(s)][1, ])
  rbind(data.table(i, type = "size", x = de$x, y = de$y), data.table(i, type = "growth", x = seq_along(ms) - 1, y = ms)) }))
sims[, wn := k$wn[i]]
env <- sims[is.finite(y), .(lo = wq(y, wn, .05), mid = sum(wn * y) / sum(wn), hi = wq(y, wn, .95)), by = .(type, x)]
in_win <- size_structure_obs$L >= from & size_structure_obs$L <= to
od <- density(size_structure_obs$L[in_win], bw = BW, from = from, to = to); obs_size <- data.table(L = od$x, density = od$y)
th <- theme_bw(base_size = 10) + theme(aspect.ratio = 1, panel.grid.minor = element_blank())
pA <- ggplot() + geom_ribbon(data = env[type == "size"], aes(x, ymin = lo, ymax = hi), fill = "grey80") +
  geom_line(data = env[type == "size"], aes(x, mid), linewidth = .8) + geom_line(data = obs_size, aes(L, density), colour = "firebrick", linewidth = .8) +
  labs(x = "Length class (cm)", y = "Density", title = expression(bold(A)~"Size distribution")) + th
amax <- 12
pB <- ggplot() + geom_jitter(data = age_structure_obs[age <= amax], aes(age, size), colour = "firebrick", size = .6, alpha = .03, width = .15, height = 0) +
  geom_ribbon(data = env[type == "growth" & x <= amax], aes(x, ymin = lo, ymax = hi), fill = "black", alpha = .25) +
  geom_line(data = env[type == "growth" & x <= amax], aes(x, mid), linewidth = .8) +
  labs(x = "Age (years)", y = "Length (mm)", title = expression(bold(B)~"Size at age")) + coord_cartesian(xlim = c(0, amax)) + scale_x_continuous(breaks = seq(0, amax, 2)) + th
of <- ogive[Sex == "F"]
pC <- ggplot(k, aes(LT)) + annotate("rect", xmin = of$L25_mm, xmax = of$L75_mm, ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = .12) +
  geom_density(aes(weight = wn), fill = "grey80", colour = "black", linewidth = .8) +
  geom_vline(xintercept = of$L50_mm, colour = "firebrick", linewidth = .8) +
  labs(x = expression("Length at maturation "*L[T]*" (mm)"), y = "Posterior density", title = expression(bold(C)~"Maturation")) + th
p <- pA + pB + pC + plot_layout(nrow = 1)
ggsave(paste0(out, ".pdf"), p, width = 17.8, height = 6.6, units = "cm", scale = 1.35)
ggsave(paste0(out, ".png"), p, width = 17.8, height = 6.6, units = "cm", scale = 1.35, dpi = 300)
## fit numbers for the caption/text
cat(sprintf("posterior L_T %.0f [%.0f, %.0f] mm vs observed female L50 %.0f (L25-L75 %.0f-%.0f); n particles %d\n",
  wq(k$LT, k$wn, .5), wq(k$LT, k$wn, .05), wq(k$LT, k$wn, .95), of$L50_mm, of$L25_mm, of$L75_mm, nrow(k)))
cat(sprintf("SRMSE size %.2f, growth %.2f (weighted medians)\n", wq(k$SRMSE_dens, k$wn, .5), wq(k$SRMSE_growth, k$wn, .5)))
