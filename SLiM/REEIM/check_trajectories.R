library(data.table)
library(ggplot2)
library(ggsci)
library(scales)

# -------------------------------
# Paths and inputs
# -------------------------------
run_root <- "./SLiM/REEIM/out/"
files <- list.files(run_root, recursive = TRUE, full.names = TRUE,
                    pattern = "H*trait_stats")

gens <- seq(9900, 10325, by = 5)

get_param <- function(key, bname) {
  as.numeric(strsplit(strsplit(bname, key)[[1]][2], "_", fixed = TRUE)[[1]][1])
}

# -------------------------------
# Load trait trajectories
# -------------------------------

#fi <- 1
traj_raw <- rbindlist(lapply(seq_along(files), function(fi) {
  fp <- files[fi]
  bname <- basename(fp)
  
  pars <- data.table(
    Ve = get_param("Ve", bname),
    sd = get_param("sd", bname),
    w = get_param("w", bname),
    f = get_param("f", bname),
    Shmax = get_param("Shmax", bname),
    SLT = get_param("SLT", bname),
    H = get_param("H", bname),
    Hx = get_param("Hx", bname),
    K = get_param("K", bname),
    theta_idx = fi
  )
  
  dt <- fread(fp)
  dt <- dt[G %in% gens, .(G, LT, g, hmax, h)]
  
  cbind(pars, dt)
}))[!is.na(G)]


# -------------------------------
# Re-center trajectories at present (G = 10125)
# -------------------------------
# Empirical reference values
ref_vals <- list(
  LT   = 169.2,
  hmax = 62,
  h    = 38,
  g    = 0.24
)

traj_centered <- traj_raw[, {
  refG <- 10125
  
  .(
    G = G,
    LT   = LT   - (mean(LT[G == refG])   - ref_vals$LT),
    g    = g    - (mean(g[G == refG])    - ref_vals$g),
    hmax = hmax - (mean(hmax[G == refG]) - ref_vals$hmax),
    h    = h    - (mean(h[G == refG])    - ref_vals$h)
  )
}, by = .(Ve, sd, K, H, Hx,theta_idx)]

traj_centered[, Year := G - 10000 + 1900]

# -------------------------------
# Apply running mean smoothing
# -------------------------------
k_smooth <- 5
traj_smoothed <- traj_centered[, .(
  G = G,
  LT   = frollmean(LT,   n = k_smooth, align = "right", na.rm = TRUE),
  g    = frollmean(g,    n = k_smooth, align = "right", na.rm = TRUE),
  hmax = frollmean(hmax, n = k_smooth, align = "right", na.rm = TRUE),
  h    = frollmean(h,    n = k_smooth, align = "right", na.rm = TRUE)
), by = .(Ve, sd, K, H, Hx, theta_idx)]

traj_smoothed[, Year := G - 10000 + 1900]


# -------------------------------
# Melt to long format
# -------------------------------
traj_long <- melt(
  traj_smoothed,
  id.vars = c("Ve", "sd", "K", "H","Hx", "theta_idx", "G", "Year"),
  measure.vars = c("h", "g", "LT")
)


traj_long[, Trait := fifelse(
  variable == "h",  "Adult growth increment (mm)",
  fifelse(variable == "LT", "Size at maturation (mm)",
          "Parental investment")
)]
# -------------------------------
# Filter baseline assumptions
# -------------------------------
traj_base <- traj_long[
  Year >= 1850 & Year <= 2225
]

# -------------------------------
# Summarise predictive envelopes
# -------------------------------
traj_sum <- traj_base[Ve==0.05 & sd==0.05 & K==10000 & Hx==1, .(
  q05 = quantile(value, 0.05, na.rm = TRUE),
  q50 = quantile(value, 0.50, na.rm = TRUE),
  q95 = quantile(value, 0.95, na.rm = TRUE)
), by = .(Year, Trait,H,theta_idx)]

# -------------------------------
# Plot trajectories with envelopes
# -------------------------------


p_pred_growth_traj <- ggplot(traj_sum, aes(x = Year,group=theta_idx,col=H)) +
  geom_line(
    aes(y = q50),
    linewidth = 0.5
  ) +
  facet_grid(Trait~., scales = "free_y") +
  geom_vline(xintercept = c(1900, 2025), linetype = 2, colour = "grey40") +
  scale_color_viridis_c(option="turbo")+
  # scale_colour_manual(
  #   values = pal_jco()(6)[c(1, 2, 4)],
  #   labels = function(x) parse(text = x),
  #   name = "Harvest rate (H)"
  # ) +
  # scale_fill_manual(
  #   values = pal_jco()(6)[c(1, 2, 4)],
  #   labels = function(x) parse(text = x),
  #   name = "Harvest rate (H)"
  # ) +
  labs(x = "Year", y = NULL) +
  theme_bw(base_size = 15) +
  theme(
    aspect.ratio = 1,
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

p_pred_growth_traj
