library(data.table)
library(ggplot2)
library(egg)

# Accepted particles from final generation
base_folder <- "./REEIM/out_REEIM/"
APOT        <- 10125     
REEIM_res   <- readRDS(file.path(base_folder, "REEIM_res.rds"))
particles   <- REEIM_res$particles
post        <- particles[[length(particles)]][keep == TRUE]

# Files corresponding to accepted particles
files_post <- post$file

# ---- Extract size distributions from simulations -----------

get_size_density <- function(file, gen, BW = 2, from = 6, to = 40) {
  
  sim_data <- fread(file, fill = Inf)
  
  size_counts <- sim_data[V1 == "countSize" & V2 == APOT, 3:ncol(sim_data)]
  size_counts <- as.numeric(size_counts)
  size_counts[is.na(size_counts)] <- 0
  
  sizes <- rep(seq_along(size_counts)-1, size_counts) ## smallest size class is 0 since size is truncated
  sim_dens <- density(sizes, bw = BW, from = from, to = to)
  
  data.table(L = sim_dens$x, density = sim_dens$y)
}

# Sample posterior simulations (for speed)
files <- paste0(base_folder,"GEN_5/","demographicData_",files_post)
pp_size <- rbindlist(lapply(files, get_size_density, gen = APOT))

# Summarise posterior predictive envelope
pp_size_sum <- pp_size[, .(
  lo = quantile(density, 0.05),
  mid = mean(density),
  hi = quantile(density, 0.95)
), by = L]

# Observed density
obs_dens <- density(size_structure_obs$L, bw = 2, from = 6, to = 40)
obs_size <- data.table(L = obs_dens$x, density = obs_dens$y)

p_size <- ggplot() +
  geom_ribbon(
    data = pp_size_sum,
    aes(L, ymin = lo, ymax = hi),
    fill = "grey80"
  ) +
  geom_line(
    data = pp_size_sum,
    aes(L, mid),
    color = "black",
    linewidth = 0.8
  ) +
  geom_line(
    data = obs_size,
    aes(L, density),
    color = "firebrick",
    linewidth = 0.8
  ) +
  theme_bw() +
  theme(aspect.ratio = 1)+
  labs(
    x = "Length class (cm)",
    y = "Density",
    title = "a) Size distribution"
  )

# ---- Extract growth curves ---------------------------------
get_growth_curve <- function(file, gen) {
  
  sim_data <- fread(file, fill = TRUE)
  
  mean_size <- sim_data[V1 == "meanSize" & V2 == gen]
  sizes     <- as.numeric(mean_size[, 3:ncol(mean_size)])
  
  data.table(
    age  = seq_along(sizes) - 1,
    size = sizes
  )
}
files <- paste0(base_folder,"GEN_5/","demographicData_",post$file)
pp_growth <- rbindlist(lapply(files, get_growth_curve, gen = APOT), idcol = "sim")

# Summarise posterior predictive growth
pp_growth_sum <- pp_growth[, .(
  lo = quantile(size, 0.05, na.rm = TRUE),
  mid = mean(size, na.rm = TRUE),
  hi = quantile(size, 0.95, na.rm = TRUE)
), by = age]

p_growth <- ggplot() +
  geom_point(
    data = age_structure_obs,
    aes(age, size),
    color = "firebrick",
    size = 2,
    alpha=0.01
  ) +
  geom_ribbon(
    data = pp_growth_sum,
    aes(age, ymin = lo, ymax = hi),
    fill = "black",
    alpha=0.25
  ) +
  geom_line(
    data = pp_growth_sum,
    aes(age, mid),
    color = "black",
    linewidth = 0.8
  ) +
  
  theme_bw() +
  theme(aspect.ratio = 1)+
  labs(
    x = "Age (years)",
    y = "Length",
    title = "b) Growth curve"
  ) +
  xlim(0,15) +
  ylim(0,450)

fig_pp <- ggarrange(p_size, p_growth, ncol = 2)
fig_pp

ggsave(
  "./Figures/PPC_REEIM.pdf",
  fig_pp,
  width = 8,
  height = 4
)

# Posterior predictive checks indicate that the inferred model reproduces the overall shape and
# scale of the observed size distribution, with the observed density closely tracking the posterior
# predictive mean across most size classes (Fig. S#). Local deviations outside the 90% pointwise
# predictive envelope reflect stochastic demographic variation and kernel density estimation error
# rather than systematic model mismatch.

# \textbf{Fig.~S\#. Posterior predictive checks for size structure and growth.}
#   Posterior predictive distributions generated from accepted ABC--SMC particles
#   (shaded regions, 90\% credible intervals) compared with observed data (red).
#   (a) Size distribution across all individuals. (b) Mean size-at-age (growth curve).
#   Solid black lines denote posterior predictive means. Agreement between observed
#   patterns and posterior predictive envelopes indicates that the inferred
#   eco-genetic model reproduces both population size structure and ontogenetic growth.
