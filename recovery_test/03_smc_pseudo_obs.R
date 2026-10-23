############################################################
##  In-silico recovery test, stronger version (PROPOSED, HEAVY):
##  full ABC-SMC on pseudo-observed datasets simulated at known H.
##
##  For each design row (recovery_test/pseudo_obs_design.csv):
##   1. simulate one pseudo-observed dataset with REEIM.slim at the
##      given H, the other parameters set to the perch posterior
##      medians (so the test is run where the real inference sits);
##   2. generation 1 reuses the existing 2,500 prior draws (no new
##      simulations), keeping the alpha closest to the pseudo-data;
##   3. generations 2..Gens propose from the weighted KDE, simulate,
##      and keep the alpha closest, exactly as run_REEIM_ABC.R does,
##      but with the pseudo-data as target.
##  The final generation gives posteriors for H, the other parameters
##  and the realised r_T, to compare with the known truth.
##
##  New simulations per dataset: 1 + (Gens - 1) * n_particles.
##  With the defaults below (Gens = 4, n_particles = 1000) that is
##  3,001 runs per dataset, ~18,000 for the 6-point design, i.e.
##  about 1.5x the 12,500 runs of the original perch ABC-SMC.
##
##  Run from the repository root (resumes where it stopped):
##    Rscript recovery_test/03_smc_pseudo_obs.R <design_id> [cores] [n_particles] [Gens]
##  Needs: slim and GNU parallel on PATH, recovery_test/output/gen1_stats.rds
############################################################

library(data.table)
library(parallel)
library(ks)
library(mvtnorm)
source("./recovery_test/R/functions.R")

args        <- commandArgs(trailingOnly = TRUE)
design_id   <- as.integer(args[1])
cores       <- if (length(args) > 1) as.integer(args[2]) else max(1, detectCores() - 1)

slim_script <- "./REEIM/REEIM.slim"
res_file    <- "./REEIM/out_REEIM/REEIM_res.rds"
smc_root    <- "./recovery_test/smc"
design_file <- "./recovery_test/pseudo_obs_design.csv"
n_particles <- if (length(args) > 2) as.integer(args[3]) else 1000
Gens        <- if (length(args) > 3) as.integer(args[4]) else 4
alpha       <- 0.2
c_kde       <- 2

prior_lim <- list(w = c(0.1, 0.6), Shmax = c(60, 120), SLT = c(100, 250),
                  f = c(0.1, 1.5), H = c(0, 0.3))
inf_names <- names(prior_lim)

## Fixed parameters: identical to the sim_parameters table in run_REEIM_ABC.R
sim_parameters <- data.table(
  w = NA, Shmax = NA, SLT = NA, f = NA, H = NA,
  K = 10000, LTLim = 250, hmaxLim = 0, gi = 0.24, h0 = 1.72, gPrior = 0.24,
  LTPrior = 169.2, hmaxPrior = 62.0, sd = 0.05, u = 1e-6, Ve = 0.05,
  gamma = 0.18, delta = 0.57, DK = 8.67, a = 3.95801e-6, b = 3.18606,
  projection = 0
)

build_slim_cmd <- function(theta, out_folder, rep = NULL) {
  sp <- copy(sim_parameters)
  sp[, (names(theta)) := as.list(theta)]
  if (!is.null(rep)) sp[, rep := rep]
  args <- paste0("-d ", names(sp), "=", as.character(unlist(sp)))
  paste("slim", paste0('-d folder=\\"', normalizePath(out_folder), '\\"'),
        paste(args, collapse = " "), slim_script)
}

run_cmds <- function(cmds) {
  f <- tempfile(fileext = ".txt")
  writeLines(cmds, f)
  system(paste("cat", f, "| parallel -j", cores), ignore.stdout = TRUE, ignore.stderr = TRUE)
}

read_folder <- function(folder) {
  files <- list.files(folder, full.names = TRUE, pattern = "^demographicData_")
  bind_draws(mclapply(files, extract_one, mc.cores = cores))
}

subset_draws <- function(x, i) list(theta = x$theta[i], dens = x$dens[i, , drop = FALSE],
                                    growth = x$growth[i, , drop = FALSE], file = x$file[i])

## ---------------------------------------------------------
## Design: H at six values, other parameters at the perch posterior
## medians (weighted, final generation of the real ABC-SMC)
## ---------------------------------------------------------
if (!file.exists(design_file)) {
  real <- readRDS(res_file)
  last <- real$particles[[length(real$particles)]][keep == TRUE]
  wmed <- function(x, w) { o <- order(x); x[o][which(cumsum(w[o]) / sum(w) >= 0.5)[1]] }
  med  <- sapply(c("w", "Shmax", "SLT", "f"), function(p) wmed(last[[p]], last$weights))
  design <- data.table(id = 1:6, H = c(0.025, 0.075, 0.125, 0.175, 0.225, 0.275),
                       t(signif(med, 3)))
  fwrite(design, design_file)
}
design <- fread(design_file)[id == design_id]
stopifnot(nrow(design) == 1)

load("./empirical_data/observed_data.RData")
gen1   <- readRDS("./recovery_test/output/gen1_stats.rds")
gc     <- growth_constants(age_structure_obs, ncol(gen1$growth))

run_dir  <- file.path(smc_root, paste0("PO_", design_id))
res_path <- file.path(run_dir, "smc_res.rds")
dir.create(file.path(run_dir, "obs"), showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------
## 1. Pseudo-observed dataset (re-run with a new rep if the burn-in
##    aborts; the abort rules can reject a draw)
## ---------------------------------------------------------
rep_i <- 1000
while (length(list.files(file.path(run_dir, "obs"), pattern = "^traitStats_")) == 0) {
  theta_obs <- sapply(unlist(design[, ..inf_names]), format, digits = 3)
  run_cmds(build_slim_cmd(theta_obs, file.path(run_dir, "obs"), rep = rep_i))
  rep_i <- rep_i + 1
  if (rep_i > 1010) stop("Pseudo-observed simulation aborted 10 times")
}
obs <- read_folder(file.path(run_dir, "obs"))
obs <- subset_draws(obs, 1)

## ---------------------------------------------------------
## 2-3. ABC-SMC against the pseudo-observed dataset
## ---------------------------------------------------------
smc <- if (file.exists(res_path)) readRDS(res_path) else
  list(truth = obs$theta, particles = list(), epsilons = numeric(), ESS = numeric(), KDE = list())

for (GEN in (length(smc$KDE) + 1):Gens) {
  if (GEN > Gens) break
  message("PO ", design_id, ": generation ", GEN)

  if (GEN == 1) {
    sims <- gen1
  } else {
    out_folder <- file.path(run_dir, paste0("GEN_", GEN))
    dir.create(out_folder, showWarnings = FALSE)
    while (length(list.files(out_folder, pattern = "^traitStats_")) < n_particles) {
      cmds <- replicate(100, {
        kd   <- smc$KDE[[GEN - 1]]
        draw <- function() setNames(as.numeric(rkde(1, kd)), kd$names)
        th   <- draw()
        while (out_of_bounds(th, prior_lim)) th <- draw()
        th   <- sapply(th, format, digits = 3)
        build_slim_cmd(th, out_folder)
      })
      run_cmds(cmds)
    }
    sims <- read_folder(out_folder)
    if (nrow(sims$theta) > n_particles) sims <- subset_draws(sims, sample(nrow(sims$theta), n_particles))
  }

  p <- copy(sims$theta)
  p[, file := sims$file]
  p[, dist := dist_to_target(obs, sims, gc)]
  smc$epsilons[GEN] <- quantile(p$dist, alpha)
  p[, keep := dist < smc$epsilons[GEN]]

  if (GEN == 1) {
    p[keep == TRUE, weights := 1 / sum(keep)]
  } else {
    prev <- smc$particles[[GEN - 1]][keep == TRUE]
    cov_mat <- smc$KDE[[GEN - 1]]$H * c_kde
    wts <- apply(as.matrix(p[keep == TRUE, ..inf_names]), 1, function(th) {
      prior_density(th, prior_lim) /
        sum(prev$weights * dmvnorm(as.matrix(prev[, ..inf_names]), mean = th, sigma = cov_mat))
    })
    p[keep == TRUE, weights := wts / sum(wts)]
  }
  smc$ESS[GEN] <- 1 / sum(p[keep == TRUE, weights]^2)
  smc$particles[[GEN]] <- p
  smc$KDE[[GEN]] <- kde(p[keep == TRUE, ..inf_names], w = p[keep == TRUE, weights * .N], binned = FALSE)
  saveRDS(smc, res_path)
}
message("PO ", design_id, " done")
