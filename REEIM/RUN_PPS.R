library(data.table)
library(parallel)
library(mvtnorm)

effective_sample_size <- function(w) 1 / sum((w / sum(w))^2)


kernel_condition_on_H <- function(particles_last_keep, HTarget, h = 0.02, min_ESS = 200) {
  stopifnot(all(c("H","weights") %in% names(particles_last_keep)))
  ker <- function(x) exp(-0.5 * (x / h)^2)
  w_star <- particles_last_keep$weights * ker(particles_last_keep$H - HTarget)
  
  if (sum(w_star) == 0) {
    # widen bandwidth aggressively if no support
    h2 <- h * 2
    w_star <- particles_last_keep$weights * exp(-0.5 * ((particles_last_keep$H - HTarget)/h2)^2)
  }
  
  w_star <- w_star / sum(w_star)
  ESS <- effective_sample_size(w_star)
  
  # adapt h upward if ESS too small
  while (ESS < min_ESS) {
    h <- h * 1.5
    w_star <- particles_last_keep$weights * exp(-0.5 * ((particles_last_keep$H - HTarget)/h)^2)
    w_star <- w_star / sum(w_star)
    ESS <- effective_sample_size(w_star)
    if (h > 0.2) break  # safety
  }
  list(weights = w_star, ESS = ESS, bandwidth = h)
}

# Resample Nθ draws by given weights (systematic resampling)
resample_index <- function(w, N) {
  w <- w / sum(w); cdf <- cumsum(w)
  u0 <- runif(1, 0, 1/N); u <- u0 + (0:(N-1))/N
  findInterval(u, c(0, cdf))  # 1-based indices
}


# Suppose 'particles' is your list from ABC; take the last gen:
base_folder  <- "./REEIM/out_REEIM/"
REEIM_res <- readRDS(file.path(base_folder, "REEIM_res.rds"))
particles   <- REEIM_res$particles
last_gen <- length(particles)
base <- particles[[last_gen]]

# keep only accepted particles with their final ABC weights
keep <- base[keep == TRUE]
#stopifnot(nrow(keep) > 0, "No retained particles found")

# optional sanity: normalize weights
keep[, weights := weights / sum(weights)]

# for no change in fishing intensity after APOT
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
  Hx          = 1,
  projection  = 1000 ## projection for ABC-SMC, so Hx also not relevant to set (only affects sims after APOT)
)

## mode, median and 95
#tab <- posterior_summary(post, params,prob = c(0.05,0.95)) From Table_REEIM_params.R
# mode, median and upper 90% quantile
H_targets <- c(0.0677, 0.127, 0.194)
Ve_grid   <- c(0.025, 0.05, 0.075)       
sd_grid   <- c(0.025, 0.05, 0.075)       
K_grid    <- c(10000, 20000,30000)       

scenario_grid <- 
  rbind(as.data.table(expand.grid(
    HTarget = H_targets,
    Ve = 0.05,
    sd = 0.05,
    K  = K_grid,
    Hx  = 1,
    stringsAsFactors = FALSE
  )),
  
  as.data.table(expand.grid(
    HTarget = H_targets,
    Ve = 0.05,
    sd = sd_grid,
    K  = 10000,
    Hx  = 1,
    stringsAsFactors = FALSE
  )),
  
  as.data.table(expand.grid(
    HTarget = H_targets,
    Ve = Ve_grid,
    sd = 0.05,
    K  = 10000,
    Hx  = 1,
    stringsAsFactors = FALSE
  )),
  as.data.table(expand.grid(
    HTarget = H_targets,
    Ve = 0.05,
    sd = 0.05,
    K  = 10000,
    Hx  = 0,
    stringsAsFactors = FALSE
  )))


scenario_grid <- scenario_grid[!duplicated(scenario_grid)]


# How many parameter draws from the conditional posterior per scenario?
Ntheta <- 120         # keep this reasonably large to propagate param uncertainty
R      <- 2           # stochastic replicates per draw (seeds)
set.seed(42)

generate_cmd <- function(theta_row, controls,rep,out_folder,
                         slim_script) {
  # theta_row is a one-row data.table with columns matching your model params,
  # e.g., w, S_hmax, S_LT, f, H, g_i, h_0, K, etc.
  
  # Build the -d defines. Start from your usual 'sim_parameters'
  sim_parameters <- list(
    K      = controls$K,
    H      = controls$HTarget,           # set scenario H explicitly
    f      = theta_row$f,
    w      = theta_row$w,
    Shmax = theta_row$Shmax,
    LTLim = theta_row$LTLim %||% 250,   # keep prior defaults if not stored
    hmaxLim = theta_row$hmaxLim %||% 0,
    gi    = theta_row$gi %||% 0.24,
    h0    = theta_row$h0 %||% 1.72,
    SLT   = theta_row$SLT,
    Ve     = controls$Ve,
    sd     = controls$sd,
    Hx     = controls$Hx,
    projection     = controls$projection,
    rep    = rep
  )
  
  # -d key=value string
  kv <- paste(paste0("-d ", names(sim_parameters), "=", sim_parameters), collapse = " ")
  
  cmd <- paste("slim ", paste0('-d folder=\\"', normalizePath(out_folder), '\\"'), kv, slim_script)
  return(cmd)
}

slim_script  <- "./REEIM/REEIM.slim"
out_root    <-  "./REEIM/scenarios/"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

gen_t0 <- 10000               # start of the decade
gen_t1 <- 10125               # end of the decade (adjust to your definition)

cores <- 10

directories <- sapply(1:nrow(scenario_grid),function(x){
  sc <- scenario_grid[x,]
  paste0("H",sc$HTarget, "_Ve",sc$Ve, "_sd",sc$sd, "_K",sc$K,"_Hx",sc$Hx)
})

# si <- 10
for(si in which(!directories %in% list.files(out_root))){
  sc <- scenario_grid[si,]
  cat(as.vector(unlist(sc)),"\n")
  
  
  # 1) Condition on H_target with a kernel; adapt bandwidth if ESS is small
  cond <- kernel_condition_on_H(keep, HTarget = sc$HTarget, h = 0.02, min_ESS = 200)
  w_star <- cond$weights
  ESS    <- cond$ESS
  h_used <- cond$bandwidth
  
  # 2) Resample Ntheta parameter draws from the conditional
  idx <- resample_index(w_star, Ntheta)
  theta_draws <- keep[idx]
  #theta_draws[,hist(H)]
  #cdf
  # 3) Run R replicates per draw
  
  out_dir <- file.path(out_root, paste0("H",sc$HTarget, "_Ve",sc$Ve, "_sd",sc$sd, "_K",sc$K,"_Hx",sc$Hx,"/"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Seeds: stable & reproducible
  
  base_seed <- 1L
  
  # i <- 1
  # rep <- 1
  cmds <- do.call(rbind,lapply(seq_len(Ntheta), function(i) {
    th <- theta_draws[i]
    
    sapply(1:R,function(rep){
      generate_cmd(
        theta_row   = th,
        controls    = list(HTarget = sc$HTarget, Ve = sc$Ve, sd = sc$sd, K = sc$K,projection=1000,Hx=sc$Hx),
        rep         = rep,
        out_folder  = out_dir,
        slim_script = slim_script
      )  
    })
    
  }))
  
  
  write.table(as.vector((cmds)),"cmds.txt",row.names = FALSE,col.names = FALSE,sep = " ",quote = FALSE)
  
  ## run the simulation in SLiM using system call
  system(paste("cat cmds.txt | parallel -j", cores),ignore.stdout = TRUE,ignore.stderr = TRUE,intern = TRUE)
}


#### === reduce fishing to 0 === ####


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
  Hx          = 0,
  projection  = 1000 
)
scenario_grid <- as.data.table(expand.grid(
  HTarget = H_targets,
  Ve = 0.05,
  sd = 0.05,
  K  = 10000,
  stringsAsFactors = FALSE
))

directories <- sapply(1:nrow(scenario_grid),function(x){
  sc <- scenario_grid[x,]
  paste0("H",sc$HTarget, "_Ve",sc$Ve, "_sd",sc$sd, "_K",sc$K,"_Hx0")
})

#si <- 1
for(si in which(!directories %in% list.files(out_root))){
  sc <- scenario_grid[si,]
  cat(as.vector(unlist(sc)),"\n")
  
  
  # 1) Condition on H_target with a kernel; adapt bandwidth if ESS is small
  cond <- kernel_condition_on_H(keep, HTarget = sc$HTarget, h = 0.02, min_ESS = 200)
  w_star <- cond$weights
  ESS    <- cond$ESS
  h_used <- cond$bandwidth
  
  # 2) Resample Ntheta parameter draws from the conditional
  idx <- resample_index(w_star, Ntheta)
  theta_draws <- keep[idx]
  #theta_draws[,hist(H)]
  #cdf
  # 3) Run R replicates per draw
  
  out_dir <- file.path(out_root, paste0("H",sc$HTarget, "_Ve",sc$Ve, "_sd",sc$sd, "_K",sc$K,"_Hx0"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Seeds: stable & reproducible
  
  base_seed <- 1L
  
  # i <- 1
  # rep <- 1
  cmds <- do.call(rbind,lapply(seq_len(Ntheta), function(i) {
    th <- theta_draws[i]
    
    sapply(1:R,function(rep){
      generate_cmd(
        theta_row   = th,
        controls    = list(HTarget = sc$HTarget, Ve = sc$Ve, sd = sc$sd, K = sc$K),
        rep         = rep,
        out_folder  = out_dir,
        slim_script = slim_script
      )  
    })
    
  }))
  
  
  write.table(as.vector((cmds)),"cmds.txt",row.names = FALSE,col.names = FALSE,sep = " ",quote = FALSE)
  
  ## run the simulation in SLiM using system call
  system(paste("cat cmds.txt | parallel -j", cores),ignore.stdout = TRUE,ignore.stderr = TRUE,intern = TRUE)
}