## Timing signature of genetic vs phenotypic change: per-decade % change in mean L_T (genetic) and in mean size at
## age 4 (phenotypic) against the decade-mean relative fishing intensity, pooled over run directories,
## with the published female Lp50 decadal rates overlaid.
## Usage: Rscript plot_timing_signature.R <out.png> <run_dir> [run_dir ...]   (Tdelta = 10000, repo root = ../../..)
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
args <- commandArgs(TRUE); out <- args[1]; dirs <- args[-1]; Tdelta <- 10000; start_year <- 1850; APOT <- Tdelta + 170
FI <- fread("../../input/plaice/fishing_intensity.csv"); setnames(FI, c("year", "FI"))
dec <- function(y) 10 * ((y - start_year) %/% 10) + start_year
rows <- list()
for (dir in dirs) {
  d <- fread(file.path(dir, "design.csv"))
  for (i in seq_len(nrow(d))) {
    tf <- Sys.glob(sprintf("%s/out/%s/traitStats_*_rep%d.txt", dir, d$set[i], d$rep[i]))
    df <- Sys.glob(sprintf("%s/out/%s/demographicData_*_rep%d.txt", dir, d$set[i], d$rep[i]))
    if (!length(tf) || !length(df)) next
    tr <- fread(tf); if (max(tr$G) < APOT) next
    tr <- tr[G >= Tdelta, .(year = start_year + G - Tdelta, LT)]
    if ("gen" %in% names(d) && d$gen[i] == "extreme") next   # h2 ~ 0.95, implausible
    if (min(tr$LT) < 50) next                                  # degenerate: L_T collapses (no fish >= 12 cm)
    dm <- fread(df, fill = Inf)[V1 == "meanSize" & V2 >= Tdelta]
    s4 <- data.table(year = start_year + dm$V2 - Tdelta, S4 = as.numeric(dm$V7))   # columns V3.. are ages 0.., so V7 = age 4
    x <- merge(tr, s4, by = "year")
    x <- x[year %% 10 == 0]                                   # decade boundaries 1850, 1860, ..., 2020
    x <- x[order(year), .(decade = head(year, -1), gen = 100 * (tail(LT, -1) / head(LT, -1) - 1),
                          phen = 100 * (tail(S4, -1) / head(S4, -1) - 1))]
    rows[[length(rows) + 1]] <- data.table(src = basename(dir), set = d$set[i], rep = d$rep[i], H = d$H[i], x)
  }
}
r <- rbindlist(rows)[FI[, .(FI = mean(FI)), by = .(decade = dec(year))], on = "decade", nomatch = 0]
r <- melt(r, measure.vars = c("gen", "phen"), variable.name = "comp", value.name = "rate")
r[, comp := factor(comp, c("gen", "phen"), c("Genetic: mean L_T", "Phenotypic: mean size at age 4"))]
## Published Lp50 (female, age 4) decadal rates, compounded over each period, at the period's mean FI
pub <- data.table(src = c("van Walraven 2010, 1905-1985", "Grift 2003, 1960-1990", "van Walraven 2010, 1985-2004"),
                  y0 = c(1905, 1960, 1985), y1 = c(1985, 1990, 2004), ratio = c(0.80, 0.875, 0.55 / 0.80))
pub[, rate := 100 * (ratio^(10 / (y1 - y0)) - 1)]
pub[, FI := mapply(function(a, b) FI[year >= a & year <= b, mean(FI)], y0, y1)]
pub <- rbindlist(lapply(levels(r$comp), function(k) copy(pub)[, comp := k]))[, comp := factor(comp, levels(r$comp))]
cs <- r[, .(rho = cor(FI, rate, method = "s")), by = comp]
cs[, lab := sprintf("Spearman rho(FI, rate) = %.2f", rho)]
fwrite(r, sub("\\.png$", ".csv", out))
p <- ggplot(r, aes(FI, rate)) + geom_hline(yintercept = 0, linetype = 3) +
  geom_point(aes(colour = src), alpha = 0.25, size = 0.8, position = position_jitter(width = 0.004, height = 0)) +
  geom_smooth(method = "loess", formula = y ~ x, colour = "black", linewidth = 0.6, se = FALSE) +
  geom_point(data = pub, aes(shape = src), size = 3, fill = "red", colour = "black") +
  scale_shape_manual(values = c(21, 24, 25)) + scale_colour_manual(values = c("#0072B2", "#E69F00", "#009E73", "#CC79A7")) +
  geom_text(data = cs, aes(x = -Inf, y = -Inf, label = lab), hjust = -0.05, vjust = -0.6, size = 3) +
  facet_wrap(~comp, nrow = 1) +
  labs(x = "Decade-mean relative fishing intensity FI", y = "Change over the decade (%)", colour = "Simulations", shape = "Published Lp50 (per decade)",
       title = "Timing signature: per-decade change vs fishing intensity (1850-2020)",
       subtitle = "Each point = one run x one decade (runs with h2 ~ 0.95 or L_T < 50 mm excluded); line = loess.\nRed symbols: published female age-4 Lp50 declines as compounded decadal rates") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", legend.box = "vertical") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
ggsave(out, p, width = 10, height = 5.5, dpi = 150); cat("Wrote", out, "\n")
print(cs); print(pub[comp == levels(r$comp)[1], .(src, rate = round(rate, 1), FI = round(FI, 2))])
print(r[, .(min_rate = round(min(rate), 1), q01 = round(quantile(rate, 0.01), 1)), by = .(comp, src)])
