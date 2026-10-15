## Virtual-ecologist analysis of survey_*.txt files written by the REEIM_vesurvey.slim copy (random females every 5 years).
## Per run: maturity-ogive L50 by survey year (logistic mature ~ length, ages >= 1), PMRN Lp50 at age 4 by decade
## (logistic first-time maturing ~ age + length on fish immature at the start of the year, pooling the decade's surveys),
## and the genetic mean L_T from traitStats. Figure: all three as % change vs ~1900, with fishing intensity and published Lp50.
## Usage: Rscript analyse_vesurvey.R <run_dir> [Tdelta=10000]   (run_dir holds design.csv with a 'label' column and out/<set>/)
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
args <- commandArgs(TRUE); setwd(args[1]); Tdelta <- if (length(args) >= 2) as.numeric(args[2]) else 10000
dir.create("fig", showWarnings = FALSE); design <- fread("design.csv")
FI <- fread("../../../input/plaice/fishing_intensity.csv"); setnames(FI, c("year", "FI"))
lp50 <- function(d, a = 4) {          # PMRN midpoint at age a from first ~ age + length
  if (sum(d$first) < 10) return(NA_real_)
  m <- glm(first ~ age + length, binomial, d); b <- coef(m); -(b[1] + b[2] * a) / b[3]
}
nboot <- 200
l50 <- function(d) { m <- glm(mature ~ length, binomial, d); b <- coef(m); -b[1] / b[2] }
res <- list()
for (i in seq_len(nrow(design))) {
  sf <- Sys.glob(sprintf("out/%s/survey_*_rep%d.txt", design$set[i], design$rep[i]))
  tf <- Sys.glob(sprintf("out/%s/traitStats_*_rep%d.txt", design$set[i], design$rep[i]))
  if (!length(sf) || !length(tf)) next
  s <- fread(sf); tr <- fread(tf)[, .(year = G - Tdelta + 1850, LT)]
  if (max(tr$year) < 2020) next
  og <- s[, .(value = l50(.SD)), by = year][, measure := "Ogive L50 (all females)"]
  s[, decade := 10 * (year %/% 10)]
  set.seed(1)   # bootstrap (fish resampled within decade) 90% interval for the PMRN Lp50
  pm <- s[mature == 0 | first == 1, {b <- replicate(nboot, lp50(.SD[sample(.N, replace = TRUE)]))
    .(value = lp50(.SD), lo = quantile(b, 0.05, na.rm = TRUE), hi = quantile(b, 0.95, na.rm = TRUE))}, by = decade][
    , .(year = decade + 5, value, lo, hi)][, measure := "PMRN Lp50 at age 4"]
  ge <- tr[year >= 1895 & year %% 5 == 0, .(year, value = LT)][, measure := "Genetic mean L_T"]
  x <- rbind(og, pm, ge, fill = TRUE)
  x[, ref := mean(value[year >= 1895 & year <= 1915], na.rm = TRUE), by = measure]
  x[, `:=`(rel = 100 * (value / ref - 1), rel_lo = 100 * (lo / ref - 1), rel_hi = 100 * (hi / ref - 1))]
  res[[length(res) + 1]] <- data.table(design[i, .(label, rep)], x)
}
res <- rbindlist(res); fwrite(res, "fig/vesurvey_series.csv")
pub <- data.table(year = c(1905, 1985, 2005, 1960, 1990), rel = c(0, -20.0, -45.0, 0, -12.5),
                  src = c(rep("van Walraven 2010 (PMRN age 4)", 3), rep("Grift 2003 (PMRN age 4), 0 at 1955-64", 2)))
lo <- max(-100, min(c(res$rel, res$rel_lo, pub$rel), na.rm = TRUE)); hi <- min(30, max(c(res$rel, res$rel_hi, 5), na.rm = TRUE)); sc <- hi - lo   # clipped: PMRN extrapolation to age 4 can explode when fish mature at ages 1-2
FI[, y := lo + FI * sc]
p <- ggplot() + geom_ribbon(data = FI, aes(year, ymin = lo, ymax = y), fill = "grey88") +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_ribbon(data = res[measure == "PMRN Lp50 at age 4"], aes(year, ymin = rel_lo, ymax = rel_hi), fill = "#D55E00", alpha = 0.2) +
  geom_line(data = res, aes(year, rel, colour = measure), linewidth = 0.7) +
  geom_point(data = res[measure == "PMRN Lp50 at age 4"], aes(year, rel, colour = measure), size = 1.2) +
  geom_line(data = pub, aes(year, rel, group = src, linetype = src)) + geom_point(data = pub, aes(year, rel, shape = src)) +
  facet_wrap(~label) + scale_colour_manual(values = c("#000000", "#0072B2", "#D55E00"), name = NULL) +
  scale_y_continuous(name = "Change since 1895-1915 (%)", sec.axis = sec_axis(~ (. - lo) / sc, name = "Relative fishing intensity (grey)")) +
  coord_cartesian(xlim = c(1895, 2022), ylim = c(lo, hi), expand = FALSE) +
  labs(x = "Year", linetype = NULL, shape = NULL, title = "Virtual ecologist: what a survey would estimate vs the genetic L_T",
       subtitle = "Random 2000 females every 5 years; ogive L50 per survey year, PMRN Lp50 (age 4) per decade (first-time maturing ~ age + length;\nshaded = 90% bootstrap interval)") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", legend.box = "vertical")
ggsave("fig/vesurvey_LT_PMRN_ogive.png", p, width = 11, height = 7.5, dpi = 150)
## windows (% change) per run and measure; PMRN uses decade midpoints (1905, 1965, 1985, 1995, 2005)
w <- function(m, a, b) res[measure == m, .(v = 100 * (value[year == b] / value[year == a] - 1)), by = .(label, rep)]$v
out <- rbindlist(lapply(unique(res$measure), function(m) {
  yrs <- if (grepl("PMRN", m)) list(c(1905, 1985), c(1965, 1995), c(1985, 2005), c(1905, 2015)) else list(c(1905, 1985), c(1960, 1990), c(1985, 2005), c(1900, 2020))
  data.table(label = res[measure == m, unique(label)], measure = m,
             w1 = w(m, yrs[[1]][1], yrs[[1]][2]), w2 = w(m, yrs[[2]][1], yrs[[2]][2]), w3 = w(m, yrs[[3]][1], yrs[[3]][2]), tot = w(m, yrs[[4]][1], yrs[[4]][2]))
}))
setnames(out, c("w1", "w2", "w3", "tot"), c("1905-85", "1960-90", "1985-2005", "total"))
fwrite(out, "fig/vesurvey_windows.csv"); print(out[order(label, measure)], digits = 3)
