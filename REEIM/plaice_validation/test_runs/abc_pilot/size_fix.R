## Run from REEIM/plaice_validation/test_runs/abc_pilot/ after copying the pilot REEIM_res.rds and GEN_1 outputs here.
## Which easy change removes the plaice size-distribution misfit? Best 20 GEN_1 particles, 2020.
suppressPackageStartupMessages(library(data.table))
load("../../data/observed_data_plaice_female.RData"); APOT <- 10170; BW <- 2
p <- fread("gen1_particles.csv")[keep == TRUE][order(dist)][1:20]
sims <- lapply(p$file, function(f) { s <- fread(file.path("GEN_1", paste0("demographicData_", f)), fill = Inf)
  cs <- as.numeric(s[V1 == "countSize" & V2 == APOT, 3:ncol(s)]); cs[is.na(cs)] <- 0; cs })
score <- function(from, to = 50, L50 = NA, k = 1) median(sapply(sims, function(cs) {
  L <- seq_along(cs) - 1; w <- if (is.na(L50)) cs else cs / (1 + exp(-(L + 0.5 - L50) / k))   # survey retention on sim only
  o <- size_structure_obs$L; o <- o[o >= from & o <= to]; od <- density(o, bw = BW, from = from, to = to)
  sl <- L[L >= from & L <= to]; sw <- w[L >= from & L <= to]
  sd_ <- density(sl, weights = sw / sum(sw), bw = BW, from = from, to = to)$y
  sqrt(mean((od$y - sd_)^2)) / sd(od$y) }))
cat("Baseline (12-50 cm):", round(score(12), 3), "\n")
for (f in c(14, 16, 18, 20)) cat("kde_from =", f, ":", round(score(f), 3), "\n")
for (L50 in c(14, 15, 16, 17, 18)) for (k in c(1, 2)) cat("retention logistic L50 =", L50, "k =", k, ":", round(score(12, L50 = L50, k = k), 3), "\n")
