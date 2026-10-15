## L_T change per grid run (1900->2020, 1955->1995) and S_LT x w heatmap.
## Usage: Rscript grid_LT_change.R <grid_dir> [Tdelta=10000]
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
args <- commandArgs(TRUE); setwd(args[1]); Tdelta <- if (length(args) >= 2) as.numeric(args[2]) else 10000
dir.create("fig", showWarnings = FALSE); design <- fread("design.csv"); APOT <- Tdelta + 171 - 1
res <- rbindlist(lapply(seq_len(nrow(design)), function(i) {
  tf <- Sys.glob(sprintf("out/%s/traitStats_*_rep%d.txt", design$set[i], design$rep[i]))
  d <- data.table(design[i], status = "aborted", LT1900 = NA_real_, LT2020 = NA_real_, LT1955 = NA_real_, LT1995 = NA_real_)
  if (!length(tf)) return(d)
  tr <- fread(tf); if (!nrow(tr) || max(tr$G) < APOT) return(d)
  tr[, year := tr$G - Tdelta + 1850]
  at <- function(y) tr[year == y, LT][1]
  d[, `:=`(status = "ok", LT1900 = at(1900), LT2020 = at(2020), LT1955 = at(1955), LT1995 = at(1995))]
  d
}))
res[, `:=`(chg_1900_2020 = 100 * (LT2020 / LT1900 - 1), chg_1955_1995 = 100 * (LT1995 / LT1955 - 1))]
fwrite(res, "fig/grid_LT_change.csv")
cell <- res[, .(n_ok = sum(status == "ok"), chg = mean(chg_1900_2020, na.rm = TRUE),
                chg5595 = mean(chg_1955_1995, na.rm = TRUE)), by = .(SLT, w)]
cell[n_ok == 0, `:=`(chg = NA, chg5595 = NA)]
cell[, lab := ifelse(n_ok == 0, "aborted", sprintf("%.1f%%%s\n(%.1f%%)", chg, ifelse(n_ok < 2, "*", ""), chg5595))]
lim <- max(abs(c(cell$chg, 12)), na.rm = TRUE)
p <- ggplot(cell, aes(factor(w), factor(SLT), fill = chg)) + geom_tile(colour = "white") +
  geom_text(aes(label = lab), size = 3) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0, limits = c(-lim, lim), na.value = "grey80") +
  labs(x = "w", y = expression(S[LT]), fill = "% change\n1900-2020",
       title = expression("Mean change in " * L[T] * " 1900" %->% "2020 (2 reps; H = 0.078, f = 0.1, " * S[hmax] * " = 117)"),
       subtitle = "Cell: mean % change 1900-2020 (in brackets 1955-1995); * = one rep aborted; grey = both aborted.\nGrift et al. 2003: ~ -12% (female age-4 Lp50, 1955-64 to ~1990)") +
  theme_bw(base_size = 10)
ggsave("fig/grid_LT_change.png", p, width = 7.5, height = 5.5, dpi = 150)
print(res[, .N, by = status])
print(cell[order(SLT, w), .(SLT, w, n_ok, chg = round(chg, 1), chg5595 = round(chg5595, 1))])
