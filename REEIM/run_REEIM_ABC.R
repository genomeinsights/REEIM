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
  cores <- args[2]
} 


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



## build SLiM command
build_slim_cmd <- function(theta, sim_parameters, script,out_folder) {
  
  keep <- names(theta)
  sim_parameters[, (keep) := as.list(theta[keep])]
  
  args <- paste0(
    "-d ", names(sim_parameters), "=",
    as.character(sim_parameters)
  )
  
  paste("slim ", paste0('-d folder=\\"', normalizePath(out_folder), '\\"'), paste(args, collapse = " "), script)
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
  
  ## Observed size density
  obs_dens <- density(size_structure_obs$L, bw = BW, from = from, to = to)
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

## =======================
##  Global configuration
## =======================

slim_script  <- "./REEIM/REEIM.slim" ## the size selectivity and historical fishing intensity are hard coded in the REEIM.slim file and come from file /R/SSF.R
base_folder  <- "./REEIM/out_REEIM/"
n_particles  <- 2500
Gens         <- 5
APOT         <- 10125        # in simulation cycles, the actual APOT is determined by the historical fishing intensity defined in REEIM.slim
alpha        <- 0.2          # cutoff for accepted particles
c            <- 2

## Observed data
load("./empirical_data/observed_data.RData")  # age_structure_obs, size_structure_obs, created by script /empirical_data/R/get_observed_data.R
LT_obs <- 169.2              # comes from litterature

## Priors

prior_lim <- list(w=c(0.1,0.6), ## cost for fecundity
                  Shmax=c(60,120), ## upper limit for hmax
                  SLT=c(100,250), ## lower limit for size at maturity
                  f=c(0.1,1.5), ## scaling factor for fecundity that determines strength of density dependence
                  H=c(0,0.3)) ## Harvest rate

## Fixed simulation parameters for ABC-SMC
## NA means the parameter is not fixed
## the REEIM.slim already has the default values and REEIM.slim also contains some default values that are not specified here
sim_parameters <- data.table(
  ## used for inference
  w           = NA,
  Shmax       = NA,
  SLT         = NA,
  f           = NA,
  H           = NA,
  ## fixed params for ABC-SMC
  K           = 10000,
  LTLim       = 250,
  hmaxLim     = 0,
  gi          = 0.24,
  h0          = 1.72,
  gPrior      = 0.24,
  LTPrior     = 169.2,
  hmaxPrior   = 62.0,
  sd          = 0.05,
  u           = 1e-6,
  Ve          = 0.05,
  gamma       = 0.18,
  delta       = 0.57,
  DK          = 8.67,
  a           = 3.95801e-6,
  b           = 3.18606,
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
      cmds <- replicate(100, {
        theta <- sapply(prior_lim, function(x) format(runif(1, x[1], x[2]),digits=3))
        build_slim_cmd(theta, sim_parameters, slim_script,out_folder)
      })
    } else {
      cmds <- replicate(100, {
        theta_cand <- rkde(1, REEIM_res$KDE[[GEN-1]])
        
        ## redraw if out of bounds     
        while(out_of_bounds(theta_cand,prior_lim)){
          theta_cand <- rkde(1, REEIM_res$KDE[[GEN-1]])
        }
        
        theta <- sapply(theta_cand, function(x) format(x,digits=3))
        names(theta) <- REEIM_res$KDE[[GEN-1]]$names
        build_slim_cmd(theta, sim_parameters, slim_script,out_folder)
      })
      
    }
    
    writeLines(cmds, "cmds.txt")
    system(paste("cat cmds.txt | parallel -j", cores),ignore.stdout = TRUE,ignore.stderr = TRUE,intern = TRUE)
    
  }
  
  
  message("Summerising generation ", GEN)
  ## Summarize
  files <- list.files(out_folder, full.names = TRUE, pattern = "demographic")
  
  REEIM_res$particles[[GEN]] <- summarise_simulations(
    files,
    par_names=par_names,
    APOT,
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
