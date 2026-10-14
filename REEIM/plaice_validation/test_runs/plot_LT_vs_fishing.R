## L_T trajectories (% change vs 1895-1905) with the relative fishing intensity on a secondary axis.
## Usage: Rscript plot_LT_vs_fishing.R <run_dir> [Tdelta=10000] [repo root=.]
## A run made with make_test_runs.py (H varies) gives fig/fig3b_LT_vs_fishing.png (one line per run, coloured by H);
## a grid made with make_grid.py (fixed H) gives fig/grid_LT_vs_fishing.png (one line per S_LT x w cell mean, faceted by w).
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
args <- commandArgs(TRUE); repo <- normalizePath(if (length(args) >= 3) args[3] else ".")
setwd(args[1]); Tdelta <- if (length(args) >= 2) as.numeric(args[2]) else 10000
start_year <- 1850; APOT <- Tdelta + 171 - 1; dir.create("fig", showWarnings = FALSE)
design <- fread("design.csv"); grid <- uniqueN(design$H) == 1
FI <- fread(file.path(repo, "REEIM/input", design$set[1], "fishing_intensity.csv"))
traj <- rbindlist(lapply(seq_len(nrow(design)), function(i) {
  tf <- Sys.glob(sprintf("out/%s/traitStats_*_rep%d.txt", design$set[i], design$rep[i]))
  if (!length(tf)) return(NULL)
  tr <- fread(tf); if (!nrow(tr) || max(tr$G) < APOT) return(NULL)
  data.table(design[i], year = start_year + tr$G - Tdelta, LT = tr$LT)
}))
traj <- traj[year >= start_year]
traj[, rel := 100 * (LT / mean(LT[year >= 1895 & year <= 1905]) - 1), by = .(set, rep)]
pm <- data.table(year = c(1905, 1985, 2004, 1960, 1990), rel = c(0, -20.0, -45.0, 0, -12.5),
                 src = c(rep("van Walraven 2010 (age 4)", 3), rep("Grift 2003 (age 4), rescaled to 0 at 1955-64", 2)))
th <- theme_bw(base_size = 10) + theme(legend.position = "bottom", legend.box = "vertical")
if (grid) {
  traj <- traj[, .(rel = mean(rel), n = .N), by = .(set, SLT, w, year)]
  lines <- geom_line(data = traj, aes(year, rel, group = SLT, colour = factor(SLT)))
  col <- scale_colour_viridis_d(name = expression(S[LT]), end = 0.9)
  fac <- facet_wrap(~w, labeller = label_both, nrow = 1)
  out <- "fig/grid_LT_vs_fishing.png"; wd <- 13
  ttl <- sprintf("Grid at H = %s: mean L_T change per S_LT x w cell (2 reps) vs relative fishing intensity", design$H[1])
} else {
  lab <- c(plaice = "Main set (g = 0.25)", plaice_gfree = "Sensitivity set (g free)")
  lines <- geom_line(data = traj, aes(year, rel, group = rep, colour = H), alpha = 0.8)
  col <- scale_colour_viridis_c(name = "H", end = 0.9)
  fac <- facet_wrap(~set, labeller = as_labeller(lab))
  out <- "fig/fig3b_LT_vs_fishing.png"; wd <- 9
  ttl <- "Simulated L_T vs relative fishing intensity and published female Lp50 declines"
}
lo <- min(c(traj$rel, pm$rel), na.rm = TRUE); hi <- max(c(traj$rel, pm$rel, 5), na.rm = TRUE); sc <- hi - lo
FI[, y := lo + intensity * sc]
p <- ggplot() + geom_ribbon(data = FI, aes(year, ymin = lo, ymax = y), fill = "grey85", alpha = 0.7) +
  geom_line(data = FI, aes(year, y), colour = "grey45", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = 3) + lines +
  geom_line(data = pm, aes(year, rel, group = src, linetype = src)) + geom_point(data = pm, aes(year, rel, shape = src)) +
  fac + col + scale_y_continuous(name = "Change in mean L_T since 1895-1905 (%)",
                                 sec.axis = sec_axis(~ (. - lo) / sc, name = "Relative fishing intensity FI(t) (grey)")) +
  coord_cartesian(ylim = c(lo, hi), xlim = c(start_year, 2020), expand = FALSE) +
  labs(x = "Year", linetype = NULL, shape = NULL, title = ttl,
       subtitle = "Realised removal(t) = H * FI(t) / 0.216, so terminal-year (2020) removal equals H") + th
ggsave(out, p, width = wd, height = 5, dpi = 150); cat("Wrote", file.path(getwd(), out), "\n")
