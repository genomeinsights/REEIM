## Run from REEIM/plaice_validation/test_runs/abc_pilot/ after copying the pilot REEIM_res.rds and GEN_1 outputs here.
## Size-distribution misfit diagnostic for the 20 best GEN_1 particles of the plaice ABC pilot.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
load("../../data/observed_data_plaice_female.RData"); Td <- 10000; APOT <- Td + 170
from <- 12; to <- 50; BW <- 2
p <- fread("gen1_particles.csv")[keep == TRUE][order(dist)][1:20]
obs <- size_structure_obs$L; od <- density(obs[obs >= from & obs <= to], bw = BW, from = from, to = to)
O <- data.table(L = od$x, obs = od$y)
dens_at <- function(sim, G) { cs <- as.numeric(sim[V1 == "countSize" & V2 == G, 3:ncol(sim)]); cs[is.na(cs)] <- 0
  s <- rep(seq_along(cs) - 1, cs); s <- s[s >= from & s <= to]; density(s, bw = BW, from = from, to = to)$y }
srmse <- function(y) sqrt(mean((O$obs - y)^2)) / sd(O$obs)
res <- list(); curves <- list(); ages <- list()
for (i in seq_len(nrow(p))) {
  sim <- fread(file.path("GEN_1", paste0("demographicData_", p$file[i])), fill = Inf)
  y2020 <- dens_at(sim, APOT); ypool <- rowMeans(sapply((APOT - 5):APOT, function(G) dens_at(sim, G)))
  res[[i]] <- data.table(i, dist = p$dist[i], SRMSE_2020 = srmse(y2020), SRMSE_pool1520 = srmse(ypool),
                         SRMSE_shift = min(sapply(seq(-5, 5, 0.5), function(s) srmse(approx(O$L + s, y2020, O$L, rule = 2)$y))),
                         best_shift = seq(-5, 5, 0.5)[which.min(sapply(seq(-5, 5, 0.5), function(s) srmse(approx(O$L + s, y2020, O$L, rule = 2)$y)))])
  curves[[i]] <- data.table(i, L = O$L, sim = y2020)
  ca <- as.numeric(sim[V1 == "countAge" & V2 == APOT, 3:ncol(sim)]); ca[is.na(ca)] <- 0
  ms <- as.numeric(sim[V1 == "meanSize" & V2 == APOT, 3:ncol(sim)])
  ages[[i]] <- data.table(i, age = seq_along(ca) - 1, n = ca / sum(ca[-1]), meanL = ms / 10)
}
res <- rbindlist(res); curves <- rbindlist(curves); ages <- rbindlist(ages)
print(res[, lapply(.SD, function(x) round(median(x), 3)), .SDcols = -1])
## where is the misfit? mean residual by length band
C <- curves[O, on = "L"][, band := cut(L, c(12, 16, 20, 25, 30, 35, 40, 50), include.lowest = TRUE)]
print(C[, .(obs_share = round(sum(obs) / sum(O$obs), 3), sim_share = round(sum(sim) / sum(O$obs) / 20, 3), sq_err_share = sum((sim - obs)^2)), by = band][, sq_err_share := round(sq_err_share / sum(sq_err_share), 2)][])
## observed age composition (age-length key applied to CPUE, i.e. age_structure_obs rows) vs simulated
oa <- age_structure_obs[, .(obs = .N), by = age][, obs := obs / sum(obs)][order(age)]
sa <- ages[age > 0, .(sim = median(n)), by = age]
print(merge(oa, sa, by = "age", all = TRUE)[age <= 15][, `:=`(obs = round(obs, 3), sim = round(sim, 3))][])
print(merge(age_structure_obs[, .(obsL = mean(size) / 10), by = age], ages[, .(simL = median(meanL)), by = age], by = "age")[age <= 12][, lapply(.SD, round, 1)])
g <- ggplot() + geom_line(data = curves, aes(L, sim, group = i), colour = "#D55E00", alpha = 0.3) +
  geom_line(data = O, aes(L, obs), linewidth = 1) + theme_bw(base_size = 10) +
  labs(x = "Length (cm)", y = "Density (12-50 cm)", title = "BTS Q3 females 2015-2020 (black) vs 20 best GEN_1 particles (2020)")
ggsave("size_diag.png", g, width = 7, height = 4, dpi = 150)
