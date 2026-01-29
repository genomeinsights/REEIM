library(ggplot2)
library(ggsci)
library(data.table)
run_root <- "./REEIM/scenarios/"
files <- list.files(run_root, recursive = TRUE, full.names = TRUE, pattern = "traitStats")

gens <- seq(9900, 10325, by = 5)

get_param <- function(key, bname) {
  as.numeric(strsplit(strsplit(bname, key)[[1]][2], "_", fixed = TRUE)[[1]][1])
}

# -------------------------------
# Load trait trajectories
# -------------------------------

fi <- 1
traj_raw <- rbindlist(lapply(seq_along(files), function(fi) {
  fp <- files[fi]
  bname <- basename(fp)
  
  pars <- data.table(
    Ve = get_param("Ve", bname),
    sd = get_param("sd", bname),
    HTarget = get_param("H", bname),
    Hx = get_param("Hx", bname),
    K = get_param("K", bname),
    thetaIdx = fi
  )
  
  dt <- fread(fp)
  dt <- dt[G %in% gens, .(G, LT, g, hmax, h)]
  
  cbind(pars, dt)
}))

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

APOT = 10125
traj_centered <- traj_raw[, {
  .(
    G = G,
    LT   = LT   - (mean(LT[G == APOT])   - ref_vals$LT),
    g    = g    - (mean(g[G == APOT])    - ref_vals$g),
    hmax = hmax - (mean(hmax[G == APOT]) - ref_vals$hmax),
    h    = h    - (mean(h[G == APOT])    - ref_vals$h)
  )
}, by = .(Ve, sd, K, HTarget, Hx,thetaIdx)]

traj_centered[, Year := G - 10000 + 1900]

#traj_centered[,table(Hx)]
# -------------------------------
# Melt to long format
# -------------------------------
traj_long <- melt(
  traj_centered,
  id.vars = c("Ve", "sd", "K", "HTarget","Hx", "thetaIdx", "G", "Year"),
  measure.vars = c("h", "g", "LT")
)


traj_long[, Trait := fifelse(
  variable == "h",  "A  Annual growth increment",
  fifelse(variable == "LT", "B  Size at maturation",
          "C  Parental investment")
)]


# -------------------------------
# Filter baseline assumptions
# -------------------------------
traj_base <- traj_long[
  Year >= 1800 & Year <= 2225
]

# -------------------------------
# Summarise predictive envelopes
# -------------------------------
traj_sum <- traj_base[Ve==0.05 & sd==0.05 & K==10000, .(
  q05 = quantile(value, 0.05, na.rm = TRUE),
  q50 = quantile(value, 0.50, na.rm = TRUE),
  q95 = quantile(value, 0.95, na.rm = TRUE)
), by = .(Year, HTarget, Hx,Trait)]

# -------------------------------
# Plot trajectories with envelopes, full figure for supplementary
# -------------------------------
traj_sum[,HTarget:=round(HTarget,3)]

traj_sum[,HT := factor(c("Low H","Typical H","High H")[match(HTarget,c(0.068,0.127,0.194))],levels=c("Low H","Typical H","High H"))]


p_pred_growth_traj <- ggplot(traj_sum, aes(x = Year, group = Hx)) +
  geom_ribbon(
    data = traj_sum[Hx==1],
    aes(ymin = q05, ymax = q95),
    alpha = 0.25, colour = NA,
    fill="firebrick"
  ) +
  geom_ribbon(
    data = traj_sum[Hx==0 & Year>=2025],
    aes(ymin = q05, ymax = q95),
    alpha = 0.25, colour = NA,
    fill="steelblue"
  ) +
  geom_line(
    data = traj_sum[Hx==1],
    aes(y = q50),
    linewidth = 1,
    color="firebrick"
  ) +
  geom_line(
    data = traj_sum[Hx==0 & Year>=2025],
    aes(y = q50),
    linewidth = 0.75,
    col="steelblue"
  ) +
  facet_grid(Trait~HT, scales = "free_y") +
  geom_vline(xintercept = c(1900, 2025), linetype = 2, colour = "grey40") +
  scale_colour_manual(
    values = pal_jco()(6)[c(1, 2, 4)],
    labels = function(x) parse(text = x),
    name = "Harvest rate (H)"
  ) +
  scale_fill_manual(
    values = pal_jco()(6)[c(1, 2, 4)],
    labels = function(x) parse(text = x),
    name = "Harvest rate (H)"
  ) +
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
ggsave(
  filename = "./Figures/pred_growth_traj.pdf",
  plot     = p_pred_growth_traj,
  width    = 10,
  height   = 8
)

# -------------------------------
# Plot trajectories with envelopes, main figure, only for median H
# -------------------------------
dt_median <- traj_sum[HT=="Typical H"]

p_pred_growth_traj_median <- ggplot(dt_median, aes(x = Year, group = Hx)) +
  geom_ribbon(
    data = dt_median[Hx==1],
    aes(ymin = q05, ymax = q95),
    alpha = 0.25, colour = NA,
    fill="firebrick"
  ) +
  geom_ribbon(
    data = dt_median[Hx==0 & Year>=2025],
    aes(ymin = q05, ymax = q95),
    alpha = 0.25, colour = NA,
    fill="steelblue"
  ) +
  geom_line(
    data = dt_median[Hx==1],
    aes(y = q50),
    linewidth = 1,
    color="firebrick"
  ) +
  geom_line(
    data = dt_median[Hx==0 & Year>=2025],
    aes(y = q50),
    linewidth = 0.75,
    col="steelblue"
  ) +
  facet_wrap(~Trait, scales = "free_y", labeller = labeller(
    Trait = c(
      "Yearly juvenile size increment (mm)" = "A  Annual growth increment",
      "Size at maturation (mm)"             = "B  Size at maturation",
      "Parental investment (unitless)"      = "C  Parental investment"
    )
  )) +
  geom_vline(xintercept = c(1900, 2025), linetype = 2, colour = "grey40") +
  scale_colour_manual(
    values = pal_jco()(6)[c(1, 2, 4)],
    labels = function(x) parse(text = x),
    name = "Harvest rate (H)"
  ) +
  scale_fill_manual(
    values = pal_jco()(6)[c(1, 2, 4)],
    labels = function(x) parse(text = x),
    name = "Harvest rate (H)"
  ) +
  labs(x = "Year", y = NULL) +
  theme_bw(base_size = 15) +
  theme(
    aspect.ratio = 1,
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
p_pred_growth_traj_median

ggsave(
  filename = "./Figures/pred_growth_traj_median.pdf",
  plot     = p_pred_growth_traj_median,
  width    = 10,
  height   = 4
)



#### ==== for SM ==== #####

traj_sum <- traj_base[Hx==1, .(
  q05 = quantile(value, 0.05, na.rm = TRUE),
  q50 = quantile(value, 0.50, na.rm = TRUE),
  q95 = quantile(value, 0.95, na.rm = TRUE)
), by = .(Year, HTarget, Hx,Trait,sd,Ve,K)]
#traj_sum[,table(sd)]
traj_sum[,sd:=factor(sd)]
traj_sum[,Ve:=factor(Ve)]
traj_sum[,K:=factor(K)]

ggplot(traj_sum[K==10000 & sd==0.05], aes(x = Year,y = q50,col=Ve,group=Ve)) +
  geom_line(
    
  ) +
  facet_grid(Trait~HTarget, scales = "free_y") +
  geom_vline(xintercept = c(1900, 2025), linetype = 2, colour = "grey40") +
  scale_colour_manual(
    values = pal_jco()(6)[c(1, 2, 4)],
    labels = function(x) parse(text = x),
    name = "Harvest rate (H)"
  ) +
  scale_fill_manual(
    values = pal_jco()(6)[c(1, 2, 4)],
    labels = function(x) parse(text = x),
    name = "Harvest rate (H)"
  ) +
  labs(x = "Year", y = NULL) +
  theme_bw(base_size = 15) +
  theme(
    aspect.ratio = 0.75,
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggplot(traj_sum[K==10000 & Ve==0.05], aes(x = Year,y = q50,col=sd,group=sd)) +
  geom_line(
    
  ) +
  facet_grid(Trait~HTarget, scales = "free_y") +
  geom_vline(xintercept = c(1900, 2025), linetype = 2, colour = "grey40") +
  scale_colour_manual(
    values = pal_jco()(6)[c(1, 2, 4)],
    labels = function(x) parse(text = x),
    name = "Harvest rate (H)"
  ) +
  scale_fill_manual(
    values = pal_jco()(6)[c(1, 2, 4)],
    labels = function(x) parse(text = x),
    name = "Harvest rate (H)"
  ) +
  labs(x = "Year", y = NULL) +
  theme_bw(base_size = 15) +
  theme(
    aspect.ratio = 0.75,
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )


ggplot(traj_sum[sd==0.05 & Ve==0.05], aes(x = Year,y = q50,col=K,group=K)) +
  geom_line(
    
  ) +
  facet_grid(Trait~HTarget, scales = "free_y") +
  geom_vline(xintercept = c(1900, 2025), linetype = 2, colour = "grey40") +
  scale_colour_manual(
    values = pal_jco()(6)[c(1, 2, 4)],
    labels = function(x) parse(text = x),
    name = "Harvest rate (H)"
  ) +
  scale_fill_manual(
    values = pal_jco()(6)[c(1, 2, 4)],
    labels = function(x) parse(text = x),
    name = "Harvest rate (H)"
  ) +
  labs(x = "Year", y = NULL) +
  theme_bw(base_size = 15) +
  theme(
    aspect.ratio = 0.75,
    strip.background = element_blank(),
    strip.text = element_text(hjust = 0),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
