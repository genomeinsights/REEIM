## Summarise and plot plaice test runs made with make_test_runs.py.
## Usage: Rscript analyse_test_runs.R <out_dir> [Tdelta=10000] [repo root=.]
## Writes <out_dir>/fig/ (four figures and test_run_summary.csv).
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
args <- commandArgs(TRUE); setwd(args[1]); Tdelta <- if (length(args) >= 2) as.numeric(args[2]) else 10000
repo <- if (length(args) >= 3) args[3] else "."
start_year <- 1850; APOT <- Tdelta + 171 - 1
load(file.path(repo, "REEIM/plaice_validation/data/observed_data_plaice_female.RData"))
LT_obs <- 219; from <- 12; to <- 50; BW <- 2
design <- fread("design.csv")
obs_L <- size_structure_obs$L[size_structure_obs$L >= from & size_structure_obs$L <= to]
od <- density(obs_L, bw = BW, from = from, to = to); obs_dens <- data.table(L = od$x, dens = od$y)
obs_age <- age_structure_obs[, .(size = mean(size), sd = sd(size), n = .N), by = age][order(age)]

res <- list(); dens <- list(); grow <- list(); traj <- list()
for (i in seq_len(nrow(design))) {
  st <- design$set[i]; rp <- design$rep[i]
  f <- Sys.glob(sprintf("out/%s/demographicData_*_rep%d.txt", st, rp))
  tf <- Sys.glob(sprintf("out/%s/traitStats_*_rep%d.txt", st, rp))
  status <- if (!length(tf)) "aborted in burn-in" else "ok"
  if (status == "ok") {
    tr <- fread(tf); if (max(tr$G) < APOT) status <- "unfinished"
  }
  if (status != "ok") { res[[i]] <- data.table(design[i], status = status); next }
  sim <- fread(f, fill = Inf)
  cs <- as.numeric(sim[V1 == "countSize" & V2 == APOT, 3:ncol(sim)]); cs[is.na(cs)] <- 0
  sizes <- rep(seq_along(cs) - 1, cs); sizes <- sizes[sizes >= from & sizes <= to]
  sd_ <- density(sizes, bw = BW, from = from, to = to)$y
  srmse_d <- sqrt(mean((obs_dens$dens - sd_)^2)) / sd(obs_dens$dens)
  ms <- as.numeric(sim[V1 == "meanSize" & V2 == APOT, 3:ncol(sim)])
  ap <- data.table(age = seq_along(ms) - 1, size_pred = ms)
  g <- merge(age_structure_obs, ap, by = "age")
  srmse_g <- g[!is.na(size_pred), sqrt(mean((size - size_pred)^2)) / sd(size)]
  ca <- as.numeric(sim[V1 == "countAge" & V2 == APOT, 3:ncol(sim)]); ca[is.na(ca)] <- 0
  q99 <- (seq_along(ca) - 1)[which(cumsum(ca) / sum(ca) >= 0.99)[1]]
  t_apot <- tr[G == APOT]
  res[[i]] <- data.table(design[i], status = status, LT = t_apot$LT, g = t_apot$g, h = t_apot$h,
                         medianSize = t_apot$medianSize, N = sum(ca), age_q99 = q99,
                         SRMSE_dens = srmse_d, SRMSE_growth = srmse_g, SAD_LT = abs(t_apot$LT - LT_obs) / LT_obs)
  dens[[i]] <- data.table(set = st, rep = rp, L = obs_dens$L, dens = sd_)
  grow[[i]] <- data.table(set = st, rep = rp, ap)
  traj[[i]] <- data.table(set = st, rep = rp, year = start_year + tr$G - Tdelta, LT = tr$LT, g = tr$g, h = tr$h)
}
res <- rbindlist(res, fill = TRUE); res[, distance := SRMSE_dens + SRMSE_growth + SAD_LT]
fwrite(res, "fig/test_run_summary.csv")
print(res[, .(n = .N), by = .(set, status)])
ok <- res[status == "ok"]
ok[, best := distance <= quantile(distance, 0.25), by = set]
dens <- rbindlist(dens)[ok[, .(set, rep, distance, best, H)], on = .(set, rep)]
grow <- rbindlist(grow)[ok[, .(set, rep, distance, best, H)], on = .(set, rep)]
traj <- rbindlist(traj)[ok[, .(set, rep, distance, best, H)], on = .(set, rep)]
lab <- c(plaice = "Main set (g = 0.25)", plaice_gfree = "Sensitivity set (g free)")
th <- theme_bw(base_size = 10)

p1 <- ggplot(dens, aes(L, dens, group = rep, colour = best)) + geom_line(alpha = 0.6) +
  geom_line(data = obs_dens, aes(L, dens), inherit.aes = FALSE, colour = "black", linewidth = 1) +
  facet_wrap(~set, labeller = as_labeller(lab)) + scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#D55E00"), labels = c("other draws", "best 25%")) +
  labs(x = "Length (cm)", y = "Density (12-50 cm)", colour = NULL, title = "Female size distribution: BTS 2015-2020 (black) vs simulated 2020") + th
ggsave("fig/fig1_size_distribution.pdf", p1, width = 9, height = 4); ggsave("fig/fig1_size_distribution.png", p1, width = 9, height = 4, dpi = 150)

p2 <- ggplot(grow[age <= 20], aes(age, size_pred, group = rep, colour = best)) + geom_line(alpha = 0.6) +
  geom_pointrange(data = obs_age[age <= 20], aes(age, size, ymin = size - sd, ymax = size + sd), inherit.aes = FALSE, size = 0.2) +
  facet_wrap(~set, labeller = as_labeller(lab)) + scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#D55E00"), labels = c("other draws", "best 25%")) +
  labs(x = "Age (years)", y = "Mean length (mm)", colour = NULL, title = "Female size at age: BTS 2015-2020 (mean ± SD) vs simulated 2020") + th
ggsave("fig/fig2_size_at_age.pdf", p2, width = 9, height = 4); ggsave("fig/fig2_size_at_age.png", p2, width = 9, height = 4, dpi = 150)

## L_T relative to its mean in 1895-1905 vs published Lp50 changes (female, age 4)
ref <- traj[year >= 1895 & year <= 1905, .(LT0 = mean(LT)), by = .(set, rep)]
traj <- traj[ref, on = .(set, rep)][, rel := 100 * (LT / LT0 - 1)]
pm <- data.table(year = c(1905, 1985, 2004, 1960, 1990), rel = c(0, -20.0, -45.0, 0, -12.5),
                 src = c("van Walraven 2010 (age 4)", "van Walraven 2010 (age 4)", "van Walraven 2010 (age 4)",
                         "Grift 2003 (age 4), rescaled to 0 at 1955-64", "Grift 2003 (age 4), rescaled to 0 at 1955-64"))
p3 <- ggplot(traj[year >= 1880], aes(year, rel, group = rep, colour = best)) + geom_line(alpha = 0.6) +
  geom_line(data = pm, aes(year, rel, group = src, linetype = src), inherit.aes = FALSE, colour = "black") +
  geom_point(data = pm, aes(year, rel, shape = src), inherit.aes = FALSE, colour = "black") +
  facet_wrap(~set, labeller = as_labeller(lab)) + scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#D55E00"), labels = c("other draws", "best 25%")) +
  labs(x = "Year", y = "Change in mean L_T since ~1900 (%)", colour = NULL, linetype = NULL, shape = NULL,
       title = "Simulated L_T vs published female Lp50 declines") + th + theme(legend.position = "bottom", legend.box = "vertical")
ggsave("fig/fig3_LT_trajectory.pdf", p3, width = 9, height = 5); ggsave("fig/fig3_LT_trajectory.png", p3, width = 9, height = 5, dpi = 150)

## Distances vs H
dl <- melt(ok, id.vars = c("set", "rep", "H"), measure.vars = c("SRMSE_dens", "SRMSE_growth", "SAD_LT", "distance"))
p4 <- ggplot(dl, aes(H, value, colour = set)) + geom_point() + facet_wrap(~variable, scales = "free_y", nrow = 1) +
  geom_vline(xintercept = 0.078, linetype = 2) + scale_colour_manual(values = c(plaice = "#0072B2", plaice_gfree = "#E69F00"), labels = lab) +
  labs(x = "Terminal harvest rate H (dashed: assessed catch/TSB 2020 = 0.078)", y = "Distance", colour = NULL) + th + theme(legend.position = "bottom")
ggsave("fig/fig4_distance_vs_H.pdf", p4, width = 10, height = 3.5); ggsave("fig/fig4_distance_vs_H.png", p4, width = 10, height = 3.5, dpi = 150)

cat("\nChange in L_T, 1900 -> 2020 (%), by set:\n")
print(traj[year == 2020, .(median = median(rel), min = min(rel), max = max(rel)), by = set])
print(traj[year == 2020 & best == TRUE, .(best_median = median(rel)), by = set])
print(ok[order(set, distance), .(set, rep, H, f, w, Shmax, SLT, LT = round(LT), g = round(g, 2), h = round(h), age_q99, SRMSE_dens = round(SRMSE_dens, 2), SRMSE_growth = round(SRMSE_growth, 2), distance = round(distance, 2))])
