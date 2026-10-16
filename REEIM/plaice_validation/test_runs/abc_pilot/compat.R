## Can the observed size distribution, growth and L_T be fitted together with the published L_T declines?
## Pilot GEN_1 (500 prior draws, H 0.06-0.10) + stress1/grid1 (wider settings), all scored at kde 16-50 cm.
## Needs abc_pilot/REEIM_res_kde16.rds + GEN_1 traitStats, and test_run_summary_kde16.csv in stress1/fig and grid1/fig
## (analyse_test_runs.R with from <- 16).
## Run from REEIM/plaice_validation/test_runs/.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
r <- readRDS("abc_pilot/REEIM_res_kde16.rds")$particles[[1]]; Td <- 10000
win <- list(`1905-85` = c(1905, 1985), `1960-90` = c(1960, 1990), `1985-2004` = c(1985, 2004), `1900-2020` = c(1900, 2020))
ch <- rbindlist(lapply(seq_len(nrow(r)), function(i) {
  f <- file.path("abc_pilot/GEN_1", paste0("traitStats_", r$file[i])); if (!file.exists(f)) return(NULL)
  tr <- fread(f)[, year := G - Td + 1850]; g <- function(y) tr[year == y, LT][1]
  data.table(r[i, .(H, SLT, w, f, Shmax, dist, SRMSE_dens, SRMSE_growth, SAD_LT)], as.data.table(lapply(win, function(v) 100 * (g(v[2]) / g(v[1]) - 1))))
}))[, src := "pilot GEN_1 (prior draws, kde 16)"]
st <- fread("stress1/fig/stress_windows_runs.csv")[status == "ok" & gen != "extreme"]
ss <- rbind(fread("stress1/fig/test_run_summary_kde16.csv")[, src2 := "stress1"], fread("grid1/fig/test_run_summary_kde16.csv")[, src2 := "grid1"], fill = TRUE)
st <- merge(st, ss[, .(src = src2, rep, SRMSE_dens, SRMSE_growth, SAD_LT)], by = c("src", "rep"))[, `:=`(dist = SRMSE_dens + SRMSE_growth + SAD_LT, src = paste(src, "(kde 16)"))]
all <- rbind(ch, st[, names(ch), with = FALSE], fill = TRUE)
fwrite(all, "abc_pilot/compat_runs.csv")
comp <- c("dist", "SRMSE_dens", "SRMSE_growth", "SAD_LT")
cat("Spearman rho between L_T change (more negative = more evolution) and fit components; negative rho = more evolution goes with WORSE fit\n")
for (s in unique(all$src)) { x <- all[src == s]; cat("\n", s, " n =", nrow(x), "\n")
  print(round(sapply(names(win), function(k) sapply(comp, function(m) cor(x[[k]], x[[m]], method = "s", use = "complete"))), 2)) }
cat("\nBest fit (min distance and components) among runs reaching at least a given decline, 1960-90 window (Grift -12.5%):\n")
for (thr in c(0, -2, -4, -6, -8, -10, -12.5)) { x <- all[`1960-90` <= thr]
  cat(sprintf("decline >= %5.1f%%: n = %3d | min dist %s | min dens %s | min growth %s | min SAD_LT %s\n", -thr, nrow(x),
      if (nrow(x)) round(min(x$dist, na.rm = TRUE), 2) else "-", if (nrow(x)) round(min(x$SRMSE_dens, na.rm = TRUE), 2) else "-",
      if (nrow(x)) round(min(x$SRMSE_growth, na.rm = TRUE), 2) else "-", if (nrow(x)) round(min(x$SAD_LT, na.rm = TRUE), 2) else "-")) }
cat("\nSame for 1905-85 (van Walraven -20%):\n")
for (thr in c(0, -5, -10, -15, -20)) { x <- all[`1905-85` <= thr]
  cat(sprintf("decline >= %4.0f%%: n = %3d | min dist %s | min SAD_LT %s | min growth %s\n", -thr, nrow(x),
      if (nrow(x)) round(min(x$dist, na.rm = TRUE), 2) else "-", if (nrow(x)) round(min(x$SAD_LT, na.rm = TRUE), 2) else "-",
      if (nrow(x)) round(min(x$SRMSE_growth, na.rm = TRUE), 2) else "-")) }
L <- melt(all, id.vars = c("src", "1960-90"), measure.vars = comp, variable.name = "component")
p <- ggplot(L, aes(`1960-90`, value, colour = src)) + geom_point(alpha = 0.5, size = 1) + geom_vline(xintercept = -12.5, linetype = 2) +
  facet_wrap(~component, scales = "free_y", nrow = 1) + theme_bw(base_size = 10) + theme(legend.position = "bottom") +
  labs(x = "Simulated change in L_T 1960-90 (%); dashed = Grift 2003", y = "Distance (lower = better fit)", colour = NULL,
       title = "Fit to present-day data vs amount of evolutionary change")
ggsave("abc_pilot/compat.png", p, width = 12, height = 4, dpi = 150)
