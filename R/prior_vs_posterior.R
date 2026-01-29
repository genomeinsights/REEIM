library(data.table)
library(ggplot2)
library(viridis)

# ---- Load ABC results --------------------------------------
base_folder <- "./REEIM/out_REEIM/"
APOT        <- 10125   
REEIM_res <- readRDS(file.path(base_folder, "REEIM_res.rds"))
particles   <- REEIM_res$particles
#post        <- particles[[length(particles)-1]][keep == TRUE]
#post[,mean(w,weights)]
#sapply(particles,function(x)mean(x$dist))
# ---- Parameter classes -------------------------------------
params_with_priors <- c("w", "f", "Shmax", "SLT", "H")
emergent_traits    <- c("g", "h", "LT")

# ---- Prior definitions (only for inferred parameters) ------
prior_lim <- list(w=c(0.1,0.6), ## cost for fecundity
                  Shmax=c(60,120), ## upper limit for hmax
                  SLT=c(100,250), ## lower limit for size at maturity
                  f=c(0.1,1.5), ## scaling factor for fecundity that determines strength of density dependence
                  H=c(0,0.3)) ## Harvest rate

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
  g  = 0.24,     # from deterministic growth calibration
  h  = 37.6,     # mean adult growth increment
  LT = 169.2     # literature-based estimate of size at maturation
  # h0 = 1.72    # fixed, discussed in text (optional)
)

# ---- Plot --------------------------------------------------

prior_post <- ggplot(df, aes(x = Value)) +
  
  
  # Posteriors
  geom_density(
    aes(weight = Weight,col=factor(Generation),group = Generation),
    
    alpha = 0.7,
    linewidth = 0.4
  ) +
  geom_density(
    data = df[ Generation==max(Generation)],
    aes(weight = Weight),
    color = "black",
    fill = "grey",
    alpha = 0.4,
    linewidth = 0.9
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
    linetype = "solid",
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
    title = "Posterior distributions and reference values used for SLiM initialization"
  )
prior_post
# ---- Save --------------------------------------------------
ggsave(
  filename = "./Figures/prior_vs_posterior.pdf",
  plot     = prior_post,
  width    = 14,
  height   = 8
)
