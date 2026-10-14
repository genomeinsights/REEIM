## Maximum-response summary for stress/grid runs: L_T % change per published window, per run and per setting (replicate mean),
## with realised h2 of L_T, ABC distance (from fig/test_run_summary.csv, so run analyse_test_runs.R first) and prior/h2 flags.
## Usage: Rscript analyse_stress.R <run_dir> [<run_dir> ...]   (Tdelta = 10000). Writes <first run_dir>/fig/stress_windows_{runs,settings}.csv
suppressPackageStartupMessages(library(data.table))
dirs <- commandArgs(TRUE); Td <- 10000
win <- list(`1905-85` = c(1905, 1985), `1960-90` = c(1960, 1990), `1985-2004` = c(1985, 2004), `1900-2020` = c(1900, 2020))
pub <- c(`1905-85` = -20.0, `1960-90` = -12.5, `1985-2004` = -31.2, `1900-2020` = NA)
runs <- rbindlist(lapply(dirs, function(dir) {
  d <- fread(file.path(dir, "design.csv")); s <- fread(file.path(dir, "fig/test_run_summary.csv"))
  if (!"gen" %in% names(d)) d[, gen := "default"]
  rbindlist(lapply(seq_len(nrow(d)), function(i) {
    base <- data.table(src = basename(dir), d[i, .(set, rep, H, f, w, Shmax, SLT, gen)])
    tf <- Sys.glob(sprintf("%s/out/%s/traitStats_*_rep%d.txt", dir, d$set[i], d$rep[i]))
    if (!length(tf)) return(base[, status := "aborted"])
    tr <- fread(tf); if (max(tr$G) < Td + 170) return(base[, status := "unfinished"])
    tr[, year := G - Td + 1850]; a <- function(y) tr[year == y, LT][1]
    ch <- sapply(win, function(v) 100 * (a(v[2]) / a(v[1]) - 1))
    w30 <- sapply(1900:1990, function(y) 100 * (a(y + 30) / a(y) - 1))
    cbind(base, status = "ok", as.data.table(as.list(ch)), steep30 = min(w30), steep30_from = (1900:1990)[which.min(w30)],
          h2_1850 = tr[G == Td, Va_LT / Vp_LT], h2_2020 = tr[G == Td + 170, Va_LT / Vp_LT],
          distance = s[rep == d$rep[i] & set == d$set[i], distance])
  }))
}), fill = TRUE)
runs[, flag := trimws(paste(ifelse(SLT > 324 | w < 0.1, "outside-prior", ""), ifelse(!is.na(h2_2020) & h2_2020 > 0.6, "h2>0.6", "")))]
ok <- runs[status == "ok"]
sett <- ok[, c(.(n = .N), lapply(.SD, mean)), by = .(src, H, SLT, w, gen),
           .SDcols = c(names(win), "steep30", "h2_1850", "h2_2020", "distance")]
sett[, flag := trimws(paste(ifelse(SLT > 324 | w < 0.1, "outside-prior", ""), ifelse(h2_2020 > 0.6, "h2>0.6", "")))]
dir.create(file.path(dirs[1], "fig"), showWarnings = FALSE)
fwrite(runs, file.path(dirs[1], "fig/stress_windows_runs.csv")); fwrite(sett, file.path(dirs[1], "fig/stress_windows_settings.csv"))
cat("Run status:\n"); print(runs[, .N, by = .(src, gen, status)])
cat("\nRealised h2 of L_T (Va_LT/Vp_LT from traitStats, at 1850 = end of burn-in and at 2020), by genetics level:\n")
print(ok[, .(n = .N, h2_1850_med = median(h2_1850), h2_1850_range = sprintf("%.2f-%.2f", min(h2_1850), max(h2_1850)),
             h2_2020_med = median(h2_2020), h2_2020_range = sprintf("%.2f-%.2f", min(h2_2020), max(h2_2020)), n_above_0.6 = sum(h2_2020 > 0.6)), by = .(src, gen)], digits = 2)
best <- function(x, k, lab) { if (!nrow(x)) return(data.table(window = k, which = lab)); r <- x[which.min(get(k))]
  data.table(window = k, which = lab, change = round(r[[k]], 1), H = r$H, SLT = r$SLT, w = r$w, gen = r$gen, h2_2020 = round(r$h2_2020, 2), distance = round(r$distance, 2), flag = r$flag) }
tab <- rbindlist(lapply(names(win), function(k) rbind(
  best(sett[distance <= 1.3], k, "setting mean, distance <= 1.3"),
  best(sett[distance <= 1.5], k, "setting mean, distance <= 1.5"),
  best(sett, k, "setting mean, any fit"),
  best(ok, k, "single run (upper bound)"), fill = TRUE)), fill = TRUE)
tab[, published := pub[window]]
cat("\nLargest L_T decline per window (%; published female age-4 Lp50 change in last column):\n"); print(tab)
cat("\nSteepest 30-year window: setting mean min =", round(min(sett$steep30), 1), "; single run min =", round(min(ok$steep30), 1),
    "starting", ok[which.min(steep30), steep30_from], "\n")
cat("\nSettings sorted by distance:\n"); print(sett[order(distance)][, lapply(.SD, function(v) if (is.numeric(v)) round(v, 2) else v)])
