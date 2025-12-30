############################################################
##  Approximate Bayesian Computation – Sequential Monte Carlo
##  Eco-genetic inference framework
##
##  Author: <Your name>
##  Project: Size-selective fishing & eco-evolutionary inference
##  Purpose:
##    - Run SLiM simulations
##    - Summarize outputs
##    - Perform ABC–SMC inference
##
##  This script is designed for reproducibility and public
##  archiving (GitHub / Zenodo / Dryad).
############################################################

# NOTE:
# Parameter values are parsed directly from SLiM output filenames.
# This ensures that ABC inference remains robust to changes in
# parameter order, naming, or future model extensions.

## =======================
##  Load required packages
## =======================

library(data.table)
library(parallel)
library(mvtnorm)
library(ks)
library(ggplot2)

## =======================
##  Helper functions
## =======================

#-----------------------------------------------------------
# Parse named parameters from SLiM output filename
#-----------------------------------------------------------
parse_slim_params <- function(filename) {
  
  fname <- basename(filename)
  fname <- sub("\\.txt$", "", fname)
  
  # Split on underscores, keep key–value pairs
  tokens <- strsplit(fname, "_")[[1]]
  
  # Keep only tokens that contain numbers
  tokens <- tokens[grepl("[0-9]", tokens)]
  
  # Split key and numeric value
  kv <- lapply(tokens, function(x) {
    key <- sub("([A-Za-z]+).*", "\\1", x)
    val <- as.numeric(sub("[^0-9\\.]", "", x))
    c(key = key, value = val)
  })
  
  kv <- do.call(rbind, kv)
  
  as.data.table(as.list(setNames(
    as.numeric(kv[, "value"]),
    kv[, "key"]
  )))
}

## build SLiM command
build_slim_cmd <- function(theta, sim_parameters, script) {
  
  sim_parameters[, c("w","S_hmax","S_LT","f","H") := as.list(theta)]
  
  args <- paste0(
    "-d ", names(sim_parameters), "=",
    as.character(sim_parameters)
  )
  
  paste("slim", paste(args, collapse = " "), script)
}

#-----------------------------------------------------------
# Summarize SLiM simulation outputs and compute ABC distance
#-----------------------------------------------------------
# Arguments:
#   files               Character vector of SLiM output files
#   gen                 Generation from which summaries are extracted
#   BW                  Bandwidth for kernel density estimation
#   from, to             Size range (cm) for density comparison
#   age_structure_obs   data.table(age, size)
#   size_structure_obs  data.table(L)
#   cores               Number of CPU cores
#
# Returns:
#   data.table with summary statistics and ABC distance
#-----------------------------------------------------------
summarise_simulations <- function(files,
                                  gen,
                                  BW = 0.5,
                                  from = 6,
                                  to = 40,
                                  age_structure_obs,
                                  size_structure_obs,
                                  cores = 1) {
  
  ## Observed size density
  obs_dens <- density(size_structure_obs$L, bw = BW, from = from, to = to)
  obs_dens <- data.table(L = obs_dens$x, obs = obs_dens$y)
  
  rbindlist(
    mclapply(files, function(file) {
      
      trait_file <- sub("demographic_data_", "trait_stats_", file)
      if (!file.exists(trait_file)) return(NULL)
      
      trait_data <- fread(trait_file)
      if (max(trait_data$G) < gen) return(NULL)
      
      # ---- Parse parameters (robust) ----
      pars <- parse_slim_params(file)
      
      # ---- Load demographic data ----
      sim_data <- fread(file)
      
      # ---- Size distribution ----
      size_counts <- sim_data[V1 == "Count_size" & V2 == gen, 3:ncol(sim_data)]
      size_counts <- as.numeric(size_counts)
      size_counts[is.na(size_counts)] <- 0
      
      sizes <- rep(seq_along(size_counts), size_counts)
      sim_dens <- density(sizes, bw = BW, from = from, to = to)$y
      
      SRMSE_dens <- sqrt(mean((obs_dens$obs - sim_dens)^2)) / sd(obs_dens$obs)
      
      # ---- Growth curve ----
      mean_size <- sim_data[V1 == "mean_size" & V2 == gen, 3:ncol(sim_data)]
      mean_size <- as.numeric(mean_size)
      
      age_pred <- data.table(
        age = seq_along(mean_size) - 1,
        size_pred = mean_size
      )
      
      growth_dt <- merge(age_structure_obs, age_pred, by = "age", all.x = TRUE)
      
      SRMSE_growth <- growth_dt[
        !is.na(size_pred),
        sqrt(mean((size - size_pred)^2)) / sd(size)
      ]
      
      # ---- Trait summaries ----
      tr <- trait_data[G == gen,
                       .(LT, g, h, hmax, median_size)]
      
      cbind(
        pars,
        tr,
        SRMSE_dens,
        SRMSE_growth,
        file = file
      )
      
    }, mc.cores = cores)
  )
}


#-----------------------------------------------------------
# Prior boundary check
#-----------------------------------------------------------
out_of_bounds <- function(theta, prior_lim) {
  any(mapply(function(x, lim) x < lim[1] | x > lim[2],
             theta, prior_lim))
}

#-----------------------------------------------------------
# Prior density
#-----------------------------------------------------------
prior_density <- function(theta, prior_lim) {
  prod(mapply(function(x, lim) dunif(x, lim[1], lim[2]),
              theta, prior_lim))
}

#-----------------------------------------------------------
# Proposal density (mixture KDE)
#-----------------------------------------------------------
proposal_density <- function(theta, centers, weights, cov_mat) {
  sum(weights * mvtnorm::dmvnorm(centers, mean = theta, sigma = cov_mat))
}

## =======================
##  Global configuration
## =======================

slim_script  <- "/SLiM/REEIM/REEIM.slim" ## the size selectivity and historical fishing intensity are hard coded in the REEIM.slim file and come from file /R/SSF.R
base_folder  <- "/SLiM/REEIM/"
n_particles  <- 2500
Gens         <- 5
APOT         <- 10125        # in simulation cycles, the actual APOT is determined by the historical fishing intensity defined in REEIM.slim
alpha        <- 0.2          # cutoff for accepted particles
cores        <- 12

## Observed data
load("./empirical_data/observed_data.RData")  # age_structure_obs, size_structure_obs, created by script /empirical_data/R/get_observed_data.R
LT_obs <- 169.2              # comes from litterature

## Priors
prior_lim <- list(
  w      = c(0.1, 0.6),
  S_hmax = c(60, 120),
  S_LT   = c(150, 300),
  f      = c(0.1, 1.5),
  H      = c(0, 0.3)
)

## Fixed simulation parameters
## NA means the parameter is not fixed
sim_parameters <- data.table(
  K           = 10000,
  H           = NA,
  f           = NA,
  w           = NA,
  S_hmax      = NA,
  LT_lim      = 250,
  hmax_lim    = 0,
  g_i         = 0.24,
  h_0         = 1.72,
  S_LT        = NA,
  g_prior     = 0.24,
  L_T_prior   = 169.2,
  h_max_prior = 62.0,
  sd          = 0.05,
  u           = 1e-6,
  Ve          = 0.05,
  gamma       = 0.18,
  delta       = 0.57,
  D_K         = 8.67,
  a           = 3.95801e-6,
  b           = 3.18606,
  print       = 0
)

## =======================
##  ABC–SMC algorithm
## =======================

particles <- vector("list", Gens)
epsilons  <- numeric(Gens)
ESS       <- numeric(Gens)
kde_prev  <- NULL

for (GEN in seq_len(Gens)) {
  
  message("Starting ABC generation ", GEN)
  
  out_folder <- file.path(base_folder, paste0("GEN_", GEN))
  dir.create(out_folder, showWarnings = FALSE)
  
  ## Run SLiM simulations
  while (length(list.files(out_folder, pattern = "trait")) < n_particles) {
    
    cmds <- replicate(100, {
      theta <- sapply(prior_lim, function(x) runif(1, x[1], x[2]))
      build_slim_cmd(theta, sim_parameters, slim_script)
    })
    
    
    writeLines(cmds, "cmds.txt")
    system(paste("cat cmds.txt | parallel -j", cores))
  }
  
  ## Summarize
  files <- list.files(out_folder, full.names = TRUE, pattern = "demographic")
  particles[[GEN]] <- summarise_simulations(
    files,
    gen_extract,
    age_structure_obs = age_structure_obs,
    size_structure_obs = size_structure_obs,
    cores = cores
  )
  
  particles[[GEN]][, SAD_LT := abs(LT - LT_obs) / LT_obs]
  particles[[GEN]][, dist := SRMSE_dens + SRMSE_growth + SAD_LT]
  
  epsilons[GEN] <- quantile(particles[[GEN]]$dist, alpha)
  particles[[GEN]][, keep := dist < epsilons[GEN]]
  
  ## Weights
  if (GEN == 1) {
    particles[[GEN]][keep == TRUE, weights := 1 / sum(keep)]
  } else {
    wts <- sapply(which(particles[[GEN]]$keep), function(i) {
      theta <- unlist(particles[[GEN]][i, .(w,S_hmax,S_LT,f,H)])
      prior_density(theta, prior_lim) /
        proposal_density(theta,
                         centers = as.matrix(particles[[GEN-1]][keep == TRUE,
                                                                .(w,S_hmax,S_LT,f,H)]),
                         weights = particles[[GEN-1]][keep == TRUE, weights],
                         cov_mat = kde_prev$H)
    })
    particles[[GEN]][keep == TRUE, weights := wts / sum(wts)]
  }
  
  ESS[GEN] <- 1 / sum(particles[[GEN]][keep == TRUE, weights]^2)
  
  ## KDE proposal
  kde_prev <- kde(
    particles[[GEN]][keep == TRUE, .(w,S_hmax,S_LT,f,H)],
    w = particles[[GEN]][keep == TRUE, weights] * n_particles,
    binned = FALSE
  )
}

## =======================
##  Save results
## =======================

saveRDS(
  list(particles = particles,
       epsilons = epsilons,
       ESS = ESS),
  file = file.path(base_folder, "ABC_results.rds")
)
