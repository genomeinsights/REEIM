## Plaice Fig. 3D: genetic mean L_T and survey PMRN Lp50 (age 4, first-time maturing known) by year, % relative
## to the 1895-1905 mean, over the PPC draws (H = 0 controls excluded). Args: ppc folder with out/, design.csv, output folder, input dir.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
a <- commandArgs(TRUE); ppc <- a[1]; dfile <- a[2]; out <- a[3]; inp <- a[4]; Td <- 10000; y0 <- 1850
d <- fread(dfile); if ("control_for" %in% names(d)) d <- d[is.na(control_for)]
lp50 <- function(x, a = 4) { if (sum(x$first) < 10) return(NA_real_); b <- coef(glm(first ~ age + length, binomial, x)); -(b[1] + b[2] * a) / b[3] }
res <- rbindlist(lapply(d$rep, function(r) {
  tf <- Sys.glob(sprintf("%s/out/traitStats_*_rep%d.txt", ppc, r)); sf <- Sys.glob(sprintf("%s/out/survey_*_rep%d.txt", ppc, r))
  if (!length(tf)) return(NULL); tr <- fread(tf)[, .(year = G - Td + y0, value = LT)][year >= 1850 & year <= 2020]
  if (max(tr$year) < 2020) return(NULL)
  tr[, rel := 100 * (value / mean(value[year >= 1895 & year <= 1905]) - 1)]
  o <- tr[, .(rep = r, measure = "genetic", year, rel)]
  if (length(sf)) { s <- fread(sf); s[, decade := 10 * (year %/% 10)]
    pm <- s[mature == 0 | first == 1, .(value = lp50(.SD)), by = decade][, .(year = decade + 5, value)]
    pm[, rel := 100 * (value / value[year == 1905] - 1)]; o <- rbind(o, pm[, .(rep = r, measure = "PMRN", year, rel)]) }
  o }))
cat("draws used:", uniqueN(res$rep), "\n")
sm <- res[, .(median = median(rel, na.rm = TRUE), q05 = quantile(rel, .05, na.rm = TRUE), q95 = quantile(rel, .95, na.rm = TRUE), n = sum(!is.na(rel))), by = .(measure, year)]
w <- dcast(sm, year ~ measure, value.var = c("median", "q05", "q95", "n"))
setnames(w, sub("^(median|q05|q95|n)_genetic$", "LT_genetic_\\1", sub("^(median|q05|q95|n)_PMRN$", "PMRN_Lp50_age4_\\1", names(w))))
fi <- fread(file.path(inp, "fishing_intensity.csv")); setnames(fi, "intensity", "fishing_intensity_rel")
w <- merge(w, fi, by = "year", all.x = TRUE)
dir.create(out, showWarnings = FALSE, recursive = TRUE); fwrite(w, file.path(out, "ppc_plaice_LT_by_year.csv"))
## panel in the perch Fig. 3 style
g <- sm[measure == "genetic"]; p <- sm[measure == "PMRN" & year <= 2015]
lo <- min(c(g$q05, p$q05, -35), na.rm = TRUE); sc <- abs(lo)
pub <- data.table(year = c(1905, 1985, 2005), rel = c(0, -20.0, -45.0), src = "van Walraven et al. 2010")
pub2 <- data.table(year = c(1960, 1990), rel = c(0, -12.5), src = "Grift et al. 2003")
pl <- ggplot() +
  geom_ribbon(data = fi[year >= 1850 & year <= 2020], aes(year, ymin = lo, ymax = fishing_intensity_rel * sc * 0.35 + lo), fill = "grey80", alpha = .6) +
  geom_ribbon(data = g, aes(year, ymin = q05, ymax = q95), fill = "firebrick", alpha = .25) +
  geom_line(data = g, aes(year, median), colour = "firebrick", linewidth = 1) +
  geom_pointrange(data = p, aes(year, median, ymin = q05, ymax = q95), colour = "grey25", size = .3) +
  geom_line(data = pub, aes(year, rel), linetype = 2, colour = "steelblue") + geom_point(data = pub, aes(year, rel), shape = 17, colour = "steelblue", size = 2.5) +
  geom_line(data = pub2, aes(year, rel), linetype = 2, colour = "darkorange") + geom_point(data = pub2, aes(year, rel), shape = 15, colour = "darkorange", size = 2.5) +
  annotate("text", 1852, lo + 0.37 * sc, label = "fishing intensity", hjust = 0, size = 3.3, colour = "grey40") +
  annotate("text", 2001, -45, label = "van Walraven\net al. 2010", hjust = 1, vjust = 0.5, size = 3.1, colour = "steelblue") +
  annotate("text", 1993, -12.5, label = "Grift et\nal. 2003", hjust = 0, size = 3.1, colour = "darkorange") +
  scale_x_continuous(limits = c(1850, 2030)) +
  geom_hline(yintercept = 0, linetype = 3) + geom_vline(xintercept = 1900, linetype = 2, colour = "grey40") +
  labs(x = "Year", y = "change relative to 1895-1905 (%)", title = expression(bold(D)~"plaice: length at maturation ("*L[T]*")")) +
  theme_bw(base_size = 15) + theme(aspect.ratio = 1, panel.grid.minor = element_blank(), plot.title = element_text(size = 13))
ggsave(file.path(out, "Fig3D_plaice_LT.pdf"), pl, width = 5.2, height = 5); ggsave(file.path(out, "Fig3D_plaice_LT.png"), pl, width = 5.2, height = 5, dpi = 180)
print(w[year %in% c(1900, 1905, 1930, 1960, 1985, 1990, 2005, 2015, 2020)], digits = 3)
