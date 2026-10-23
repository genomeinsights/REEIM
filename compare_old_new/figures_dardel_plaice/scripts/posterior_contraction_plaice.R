library(data.table); library(ggplot2); library(gridExtra)
library(data.table)
library(ggplot2)
library(viridis)

# ---- Load ABC results --------------------------------------
base_folder <- "/private/tmp/claude-539526166/-Users-petrikem-gitlab-REEIM/f7799891-e671-4825-bad9-a218acb19f48/scratchpad/final/res_plaice/"
APOT        <- 10125   
REEIM_res <- readRDS(file.path(base_folder, "REEIM_res.rds"))
particles   <- REEIM_res$particles
post <- particles[[length(particles)]][keep == TRUE]
#post        <- particles[[length(particles)-1]][keep == TRUE]
#post[,mean(w,weights)]
#sapply(particles,function(x)mean(x$dist))
# ---- Parameter classes -------------------------------------
params_with_priors <- c("w", "f", "Shmax", "SLT", "H")
emergent_traits    <- c("g", "h", "LT")

# ---- Prior definitions (only for inferred parameters) ------
## read from the priors used by the ABC run, so the plotted priors match it
priors_csv <- fread(file.path(Sys.getenv("REEIM_INPUT", "./REEIM/input/perch"), "priors.csv"))
prior_lim <- setNames(lapply(seq_len(nrow(priors_csv)), function(i) c(priors_csv$lower[i], priors_csv$upper[i])),
                      priors_csv$parameter)

# ---- Expression labels -------------------------------------
param_labels <- c(
  w      = "w",
  f      = "f",
  Shmax = "S[h[max]]",
  SLT   = "S[L[T]]",
  H      = "H",
  g      = "g",
  h      = "h",
  LT     = "L[T]"
)

# ---- Sample from priors ------------------------------------
sample_prior <- function(name, n = 50000) {
  lim <- prior_lim[[name]]
  data.table(
    Generation = 0,
    Parameter = name,
    Value     = runif(n, lim[1], lim[2]),
    Weight    = 1,
    Source    = "Posterior"
  )
}

# ---- Build plotting data -----------------------------------
df_list <- list()

# Parameters with priors: prior + posterior
for (p in params_with_priors) {
  df_list[[length(df_list) + 1]] <- sample_prior(p)
  
  for(G in 1:length(particles)){
    df_list[[length(df_list) + 1]] <- data.table(
      Generation=G,
      Parameter = p,
      Value     = particles[[G]][which(keep)][[p]],
      Weight    = particles[[G]][which(keep),weights],
      Source    = "Posterior"
    )  
  }
}

# Emergent traits: posterior only
p <- "g"
#G <- 1
for (p in emergent_traits) {
  for(G in 1:length(particles)){
    df_list[[length(df_list) + 1]] <- data.table(
      Generation=G,
      Parameter = p,
      Value     = particles[[G]][which(keep)][[p]],
      Weight    = particles[[G]][which(keep),weights],
      Source    = "Posterior"
    )  
  }
}

facet_levels <- c(
  "H",
  "S[h[max]]",
  "S[L[T]]",
  "f",
  "w",
  "L[T]",
  "g",
  "h"
)


df <- rbindlist(df_list,fill = TRUE)

df[, Par := factor(
  param_labels[Parameter],
  levels = facet_levels
)]

ref_values <- list(
  g  = 0.25,     # plaice gPrior
  h  = 50.0,     # plaice fitted h
  LT = 219.0     # plaice LT_obs
  # h0 = 1.72    # fixed, discussed in text (optional)
)

# ---- Plot --------------------------------------------------

prior_post <- ggplot(df, aes(x = Value)) +
  
  
  # Posteriors
  # geom_density(
  #   data = df[ !Generation %in% c(0,5)],
  #   aes(weight = Weight,col=factor(Generation),group = Generation),
  #   alpha = 0.7,
  #   linewidth = 0.4
  # ) +
  geom_density(
    data = df[ Generation==max(Generation)],
    aes(weight = Weight),
    color = "black",
    fill = "grey",
    alpha = 0.4,
    linewidth = 0.9
  ) +
  geom_density(
    data = df[Generation == 0],
    aes(weight = Weight),
    color = "black",
    alpha = 0.4,
    linewidth = 0.9,
    linetype=2
  )+
  scale_color_viridis_d(
    option = "turbo",
    direction = -1,
    name = "ABC \ngeneration"
  )+
  # Reference lines for emergent traits
  
  geom_vline(
    data = data.table(
      Parameter = names(ref_values),
      Value     = unlist(ref_values),
      Par       = factor(
        param_labels[names(ref_values)],
        levels = facet_levels
      )
    ),
    aes(xintercept = Value),
    linetype = 2,
    color = "firebrick",
    linewidth = 0.8
  ) +
  
  facet_wrap(
    ~ Par,
    scales   = "free",
    ncol     = 4,
    labeller = label_parsed
  ) +
  
  theme_bw(base_size = 15) +
  theme(
    aspect.ratio      = 1,
    panel.grid.minor  = element_blank(),
    strip.background = element_blank()
  ) +
  labs(
    x = "Parameter value",
    y = "Density",
    title = "Prior versus posterior distributions"
  )
prior_post
# ---- Save --------------------------------------------------
ggsave(
  filename = "./compare_old_new/figures_dardel_plaice/prior_vs_posterior_plaice.pdf",
  plot     = prior_post,
  width    = 14,
  height   = 8
)


df[Parameter=="H" & Generation %in% c(0,5)]


HRate <- fread("./empirical_data/harvestrate.csv")

HRate[,Value:=harvestperc/mean(harvestperc)*df[Parameter=="H" & Generation %in% c(5),mean(Value)]]

p1 <- ggplot(df[Parameter=="H" & Generation %in% c(0,5)], aes(x = Value*100)) +
  
  geom_density(
    data=df[Parameter=="H" & Generation %in% 5],
    aes(weight = Weight,group = Generation),
    fill="grey",
    alpha = 0.7,
    linewidth = 0.6
  )+
  # Posteriors
  geom_density(
    data=df[Parameter=="H" & Generation %in% c(0)],
    aes(weight = Weight,group = Generation),
    alpha = 0.7,
    linewidth = 0.6,
    linetype=2
  )+
 

  theme_bw(base_size = 15) +
  theme(
    aspect.ratio      = 1,
    panel.grid.minor  = element_blank(),
    strip.background = element_blank()
  ) +
  labs(
    x = "Harvest rate (% removed annually)",
    y = "Density",
    title = "A Prior and posterior distributions \nof harvest intensity (H)"
  )
scale01 <- function(x) (x - mean(x)) / sd(x)

p2 <- ggplot() +
  geom_density(aes(x = scale01(post$H)), fill="grey", alpha=0.6) +
  geom_density(aes(x = scale01(HRate$harvestperc)), linetype=2,linewidth = 0.6) +
  labs(x="Standardized harvest intensity", y="Density",title =  "B Qualitative agreement between \ninferred and empirical harvest intensity")+
  theme_bw(base_size = 15) +
  theme(
    aspect.ratio      = 1,
    panel.grid.minor  = element_blank(),
    strip.background = element_blank()
  ) 

fig_posteriorH <- grid.arrange(p1,p2,nrow=1)

ggsave(
  filename = "./compare_old_new/figures_dardel_plaice/posterior_H_plaice.pdf",
  plot     = fig_posteriorH,
  width    = 12,
  height   = 6
)

#### for heritabilities and CVs ####
base_folder <- "/private/tmp/claude-539526166/-Users-petrikem-gitlab-REEIM/f7799891-e671-4825-bad9-a218acb19f48/scratchpad/final/res_plaice/"
APOT        <- 10125     
REEIM_res   <- readRDS(file.path(base_folder, "REEIM_res.rds"))
particles   <- REEIM_res$particles
post <- particles[[length(particles)]][keep == TRUE]
post        <- particles[[length(particles)]][keep == TRUE]

files <- paste0(base_folder,"GEN_5/","traitStats_",post$file)

dt <- fread(files[1])

summarize_century <- function(dt, t0, t1) {
  stopifnot(is.data.table(dt), "G" %in% names(dt))
  x <- dt[G >= t0 & G <= t1]
  if (nrow(x) == 0L) stop("No rows in the requested window.")
  
  # Per-generation derived metrics
  x[, `:=`(
    h2_LT   = Va_LT   / Vp_LT,
    IA_LT   = Va_LT   / (LT^2),
    h2_g    = Va_g    / Vp_g,
    IA_g    = Va_g    / (g^2),
    h2_hmax = Va_hmax / Vp_hmax,
    IA_hmax = Va_hmax / (hmax^2)
  )]
  
  # Helper: log-change between endpoints
  get_r <- function(v) {
    v0 <- x[G == t0, get(v)][1]
    v1 <- x[G == t1, get(v)][1]
    if (is.na(v0) || is.na(v1) || v0 <= 0 || v1 <= 0) return(NA_real_)
    log(v1) - log(v0)
  }
  
  Tlen <- length(unique(x$G)) - 1L  # should be 100 if generations are contiguous
  
  # Log and percent changes
  r_LT    <- get_r("LT")
  r_g     <- get_r("g")
  r_hmax  <- get_r("hmax")
  
  pct <- function(r) 100 * (exp(r) - 1)
  
  out <- data.table(
    trait     = c("LT","g","hmax"),
    r_total   = c(r_LT, r_g, r_hmax),
    r_per_gen = c(r_LT, r_g, r_hmax) / Tlen,
    pct_total = pct(c(r_LT, r_g, r_hmax)),
    pct_per_gen = pct(c(r_LT, r_g, r_hmax) / Tlen),
    h2_bar    = c(mean(x$h2_LT,   na.rm=TRUE)*100,
                  mean(x$h2_g,    na.rm=TRUE)*100,
                  mean(x$h2_hmax, na.rm=TRUE)*100),
    IA_bar    = c(mean(x$IA_LT,   na.rm=TRUE),
                  mean(x$IA_g,    na.rm=TRUE),
                  mean(x$IA_hmax, na.rm=TRUE)),
    VA_bar    = c(mean(x$Va_LT,   na.rm=TRUE),
                  mean(x$Va_g,    na.rm=TRUE),
                  mean(x$Va_hmax, na.rm=TRUE))
  )
  
  # Also return the per-generation table with derived quantities for plotting
  list(
    per_generation = x,
    window_summary = out,
    t0 = t0, t1 = t1, Tlen = Tlen
  )
}
fp <- files[1]
r <- 1
run_summaries <- rbindlist(lapply(1:nrow(post),function(r){
  
  fp <- paste0(base_folder,"GEN_5/","traitStats_",post$file[r])
  cat(fp,"\n")
  bname <<- basename(fp)
  
  pars <- data.table(
    Shmax = get_param("Shmax",bname),
    SLT   = get_param("SLT",bname),
    w      = get_param("w",bname),
    hlim  = get_param("hlim",bname),
    LTlim = get_param("LTlim",bname),
    f      = get_param("f",bname),
    Ve    = get_param("Ve",bname),
    sd    = get_param("sd",bname),
    H      = get_param("H",bname),
    K      = get_param("K",bname)
  )
  
  dt <- fread(fp)
  if(max(dt$G)==10126){
    summary <- summarize_century(dt,t0,t1)
    cbind(pars,weight=post$weights[r],summary$window_summary)  
  }else{
    return(NULL)
  }
}))

#dt <- melt(run_summaries,measure.vars = c("pct_total","h2_bar","IA_bar"))

#run_summaries[,plot(pct_per_gen,IA_bar)]

prior_post <- ggplot(dt, aes(y = abs(value),col=trait)) +
  
  
  # Posteriors
  # geom_density(
  #   data = df[ !Generation %in% c(0,5)],
  #   aes(weight = Weight,col=factor(Generation),group = Generation),
  #   alpha = 0.7,
  #   linewidth = 0.4
  # ) +
  geom_boxplot(
    aes(weight = weight),
    alpha = 0.4,
    linewidth = 0.9
  ) +
  facet_wrap(
   . ~ variable,
    scales   = "free",
    ncol     = 3,
    labeller = label_parsed
  ) +
  
  theme_bw(base_size = 15) +
  theme(
    aspect.ratio      = 1,
    panel.grid.minor  = element_blank(),
    strip.background = element_blank()
  ) +
  labs(
    x = "Parameter value",
    y = "Density",
    title = "Prior versus posterior distributions"
  )
prior_post
