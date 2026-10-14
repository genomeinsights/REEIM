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
##  
## You can add working directory using Rscript from the commandline, 
## otherwise the folder from which Rscript is run from will default as wd.
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

args <- commandArgs(trailingOnly = TRUE) 
if(length(args)>0){
  setwd(args[1])
  cores <- as.integer(args[2])
} 
if (!exists("cores")) cores <- max(1, parallel::detectCores() - 1)

## Folder with the input files describing the study system (see ./REEIM/input/README.md).
## Pass it as the third command-line argument to run another system.
input_dir <- if (length(args) > 2) args[3] else "./REEIM/input/perch"

## Output folder (fourth argument), so that several systems can run from one copy of the repository
base_folder <- if (length(args) > 3) args[4] else "./REEIM/out_REEIM/"


#setwd("~/REEIM/REEIM/")
## =======================
##  Helper functions
## =======================

#-----------------------------------------------------------
# Parse named parameters from SLiM output filename
#-----------------------------------------------------------


get_param <- function(key, bname) {
  as.numeric(strsplit(strsplit(bname, key)[[1]][2], "_", fixed = TRUE)[[1]][1])
}



#-----------------------------------------------------------
# Read a name,value input file into a named list
#-----------------------------------------------------------
read_settings <- function(file) {
  dt <- read.csv(file, colClasses = "character", comment.char = "#", strip.white = TRUE)
  setNames(as.list(dt$value), dt$name)
}

## build SLiM command
build_slim_cmd <- function(theta, sim_parameters, script,out_folder) {
  
  keep <- names(theta)
  sim_parameters[, (keep) := as.list(theta[keep])]
  
  args <- paste0(
    "-d ", names(sim_parameters), "=",
    as.character(sim_parameters)
  )
  
  paste("slim ",
        paste0('-d folder=\\"', normalizePath(out_folder), '\\"'),
        paste0('-d inputDir=\\"', normalizePath(input_dir), '\\"'),
        paste(args, collapse = " "), script)
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

#file <- files[1]


#file <- "./SLiM/REEIM/out//GEN_1/demographicData_K10000_Ve0.05_sd0.05_Shmax126_SLT163_w0.461_f0.568_H0.187_Hx1_rep1.txt"
#file <- paste0(out_folder,"/demographicData_",particles[[1]][which.min(SRMSE_dens),file])

summarise_simulations <- function(files,
                                  par_names = c("H","f","w","Shmax","SLT"),
                                  APOT=10125,
                                  BW = 2,
                                  from = 6,
                                  to = 40,
                                  age_structure_obs,
                                  size_structure_obs) {
  
  ## Observed size density, normalised over the comparison window only, so that
  ## size classes the survey gear under-samples don't shift the curve (same for simulated sizes below)
  obs_L <- size_structure_obs$L[size_structure_obs$L >= from & size_structure_obs$L <= to]
  obs_dens <- density(obs_L, bw = BW, from = from, to = to)
  obs_dens <- data.table(L = obs_dens$x, obs = obs_dens$y)
  
  rbindlist(
    lapply(files, function(file) {
      
      trait_file <- sub("demographicData_", "traitStats_", file)
      if (!file.exists(trait_file)) return(NULL)
      
      trait_data <- fread(trait_file)
      if (max(trait_data$G) < APOT) return(NULL)
      
      # ---- Parse parameters (robust) ----
      pars <- sapply(par_names,get_param,file)
      
      # ---- Load demographic data ----
      sim_data <- fread(file,fill = Inf)
      
      # ---- Size distribution ----
      size_counts <- sim_data[V1 == "countSize" & V2 == APOT, 3:ncol(sim_data)]
      size_counts <- as.numeric(size_counts)
      size_counts[is.na(size_counts)] <- 0
      
      sizes <- rep(seq_along(size_counts)-1, size_counts) ## smallest size class is 0 since size is truncated
      sizes <- sizes[sizes >= from & sizes <= to]
      sim_dens <- density(sizes, bw = BW, from = from, to = to)$y
      
      SRMSE_dens <- sqrt(mean((obs_dens$obs - sim_dens)^2)) / sd(obs_dens$obs)
    
      
      # ---- Growth curve ----
      mean_size <- sim_data[V1 == "meanSize" & V2 == APOT, 3:ncol(sim_data)]
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
      tr <- trait_data[G == APOT,.(LT, g, h, hmax, medianSize)]
      
      data.table(
        t(pars),
        tr,
        SRMSE_dens,
        SRMSE_growth,
        file = gsub("demographicData_","",basename(file))
      )
      
    })
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

#-----------------------------------------------------------
# Draw one candidate from the same mixture proposal_density() evaluates:
# pick a previous particle by weight, add MVN(0, cov_mat) noise
#-----------------------------------------------------------
draw_proposal <- function(centers, weights, cov_mat) {
  j <- sample.int(nrow(centers), 1, prob = weights)
  setNames(centers[j, ] + as.vector(mvtnorm::rmvnorm(1, sigma = cov_mat)), colnames(centers))
}

## =======================
##  Global configuration
## =======================

slim_script  <- "./REEIM/REEIM.slim" ## the size selectivity and historical fishing intensity are read from input_dir (created by /R/SSF.R for perch)
n_particles  <- 2500
n_batch      <- max(100, 2 * cores)  # simulations started per round; keeps all cores busy on large nodes
Gens         <- 5

## System description (input_dir)
model_parameters <- read_settings(file.path(input_dir, "parameters.csv"))   # constants read by REEIM.slim
abc_settings     <- read_settings(file.path(input_dir, "abc_settings.csv"))
fishing          <- fread(file.path(input_dir, "fishing_intensity.csv"))
priors           <- fread(file.path(input_dir, "priors.csv"))

## Cycle at which the observed (present-day) data are compared with the simulations:
## the last year of the fishing trajectory, i.e. burn-in + number of fishing years - 1
## (10125 for perch: 10000 burn-in cycles and 126 rows from 1900 to 2025)
APOT         <- as.numeric(model_parameters$Tdelta) + nrow(fishing) - 1
alpha        <- 0.2          # cutoff for accepted particles
c            <- 2

## Observed data
load(abc_settings$observed_data)  # age_structure_obs, size_structure_obs; for perch created by script /empirical_data/R/get_observed_data.R
LT_obs <- as.numeric(abc_settings$LT_obs)              # comes from litterature
kde_bw   <- as.numeric(abc_settings$kde_bw)            # bandwidth and range (cm) of the size-distribution comparison
kde_from <- as.numeric(abc_settings$kde_from)
kde_to   <- as.numeric(abc_settings$kde_to)

## Priors (uniform; the order of the rows in priors.csv sets the parameter order)

prior_lim <- setNames(lapply(seq_len(nrow(priors)), function(i) c(priors$lower[i], priors$upper[i])),
                      priors$parameter)

## Simulation parameters passed to SLiM for ABC-SMC
## NA means the parameter is inferred and set from theta
## All fixed parameters (K, growth, genetics, length-weight, burn-in limits, ...) are read
## by REEIM.slim from parameters.csv in input_dir; add a column here only to override one
sim_parameters <- data.table(
  ## used for inference
  t(setNames(rep(NA, nrow(priors)), priors$parameter)),
  projection  = 0 ## no projection for ABC-SMC, so Hx also not relevant to set (only affects sims after APOT)
)

## =======================
##  ABC–SMC algorithm
## =======================


par_names <- names(prior_lim)

#GEN = 2

## read in existing data 
if(file.exists(file.path(base_folder, "REEIM_res.rds"))){
  REEIM_res <- readRDS(file.path(base_folder, "REEIM_res.rds"))
  
} else {
  message("No previous REEIM file exists, generating new")
  REEIM_res <- list(particles = NULL,
                    epsilons = NULL,
                    ESS = NULL,
                    KDE = NULL)
}

#GEN = 2
no_gens_done <- length(REEIM_res$KDE)
if(no_gens_done==Gens) stop("Final generation has already been done")

for (GEN in (length(REEIM_res$KDE)+1):Gens) {
      
  message("Starting ABC generation ", GEN)
  out_folder <- file.path(base_folder, paste0("GEN_", GEN))

  message("Output folder: ",normalizePath(out_folder))
  dir.create(out_folder, showWarnings = FALSE,recursive = TRUE)
  
  ## Run SLiM simulations
  while (length(list.files(out_folder, pattern = "trait")) < (n_particles)) { 
    
    message(length(list.files(out_folder, pattern = "trait"))," out of ", n_particles, " particles are done")
    
    if(GEN==1){
      cmds <- replicate(n_batch, {
        theta <- sapply(prior_lim, function(x) format(runif(1, x[1], x[2]),digits=3))
        build_slim_cmd(theta, sim_parameters, slim_script,out_folder)
      })
    } else {
      prev_centers <- as.matrix(REEIM_res$particles[[GEN-1]][keep == TRUE, ..par_names])
      prev_weights <- REEIM_res$particles[[GEN-1]][keep == TRUE, weights]
      prev_cov     <- REEIM_res$KDE[[GEN-1]]$H*c   # same kernel as the importance weights below
      cmds <- replicate(n_batch, {
        theta_cand <- draw_proposal(prev_centers, prev_weights, prev_cov)
        
        ## redraw if out of bounds     
        while(out_of_bounds(theta_cand,prior_lim)){
          theta_cand <- draw_proposal(prev_centers, prev_weights, prev_cov)
        }
        
        theta <- sapply(theta_cand, function(x) format(x,digits=3))
        names(theta) <- par_names
        build_slim_cmd(theta, sim_parameters, slim_script,out_folder)
      })
      
    }
    
    cmd_file <- file.path(out_folder, "cmds.txt")  # per run, so parallel ABCs don't overwrite each other's list
    writeLines(cmds, cmd_file)
    system(paste("cat", shQuote(cmd_file), "| parallel -j", cores),ignore.stdout = TRUE,ignore.stderr = TRUE,intern = TRUE)
    
  }
  
  
  message("Summerising generation ", GEN)
  ## Summarize
  files <- list.files(out_folder, full.names = TRUE, pattern = "demographic")
  
  REEIM_res$particles[[GEN]] <- summarise_simulations(
    files,
    par_names=par_names,
    APOT,
    BW = kde_bw,
    from = kde_from,
    to = kde_to,
    age_structure_obs = age_structure_obs,
    size_structure_obs = size_structure_obs
  )
  
  ## keep 
  if(REEIM_res$particles[[GEN]][,.N]>n_particles) REEIM_res$particles[[GEN]] <- REEIM_res$particles[[GEN]][sample(1:.N,n_particles,replace = FALSE)] # keep only n_particles
  
  
  REEIM_res$particles[[GEN]][, SAD_LT := abs(LT - LT_obs) / LT_obs]
  REEIM_res$particles[[GEN]][, dist := SRMSE_dens + SRMSE_growth + SAD_LT]
  
  REEIM_res$epsilons[GEN] <- quantile(REEIM_res$particles[[GEN]]$dist, alpha)
  REEIM_res$particles[[GEN]][, keep := dist < REEIM_res$epsilons[GEN]]
  
  
  ## Weights
  #i <- 1
  if (GEN == 1) {
    REEIM_res$particles[[GEN]][keep == TRUE, weights := 1 / sum(keep)]
  } else {
    wts <- sapply(which(REEIM_res$particles[[GEN]]$keep), function(i) {
      theta <- unlist(REEIM_res$particles[[GEN]][i, ..par_names])
      prior_density(theta, prior_lim) /
        proposal_density(theta,
                         centers = as.matrix(REEIM_res$particles[[GEN-1]][keep == TRUE, ..par_names]),
                         weights = REEIM_res$particles[[GEN-1]][keep == TRUE, weights],
                         cov_mat = REEIM_res$KDE[[GEN-1]]$H*c)
    })
    REEIM_res$particles[[GEN]][keep == TRUE, weights := wts / sum(wts)]
  }
  
  REEIM_res$ESS[GEN] <- 1 / sum(REEIM_res$particles[[GEN]][keep == TRUE, weights]^2)
  #names(REEIM_res$ESS)

  ## KDE proposal
  message("Kernel density estimation ")
  
  ## This takes some time, so if the kde already exist, read from file
  REEIM_res$KDE[[GEN]] <- kde(
    REEIM_res$particles[[GEN]][keep == TRUE, ..par_names],
    w = REEIM_res$particles[[GEN]][keep == TRUE, weights*.N],
    binned = FALSE)
  
  
  ##  Save results
  saveRDS(REEIM_res,file = file.path(base_folder, "REEIM_res.rds"))
  
}
