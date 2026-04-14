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
  
  # Assumes these columns exist in traitStats:
  # LT, g, hmax, h, Vp_LT, Vp_g, Vp_hmax, Vp_h
  dt <- dt[G %in% gens, .(G, LT, g, hmax, h, Vp_LT, Vp_g, Vp_hmax)]
  
  cbind(pars, dt)
}))

# -------------------------------
# Re-center trajectories at present (G = 10125)
# -------------------------------
ref_vals <- list(
  LT   = 169.2,
  hmax = 62,
  h    = 38,
  g    = 0.24
)

APOT <- 10125

traj_centered <- traj_raw[, {
  .(
    G = G,
    LT   = LT   - (mean(LT[G == APOT])   - ref_vals$LT),
    g    = g    - (mean(g[G == APOT])    - ref_vals$g),
    hmax = hmax - (mean(hmax[G == APOT]) - ref_vals$hmax),
    Vp_LT   = Vp_LT,
    Vp_g    = Vp_g,
    Vp_hmax = Vp_hmax
  )
}, by = .(Ve, sd, K, HTarget, Hx, thetaIdx)]

traj_centered[, Year := G - 10000 + 1900]

# -------------------------------
# Existing long-format code
# -------------------------------
traj_long <- melt(
  traj_centered,
  id.vars = c("Ve", "sd", "K", "HTarget", "Hx", "thetaIdx", "G", "Year"),
  measure.vars = c("hmax","g", "LT")
)

traj_long[, Trait := fifelse(
  variable == "h",  "A  Annual growth increment",
  fifelse(variable == "LT", "B  Size at maturation",
          "C  Parental investment")
)]

traj_base <- traj_long[Year >= 1800 & Year <= 2225]

traj_sum <- traj_base[Ve == 0.05 & sd == 0.05 & K == 10000, .(
  q05 = quantile(value, 0.05, na.rm = TRUE),
  q50 = quantile(value, 0.50, na.rm = TRUE),
  q95 = quantile(value, 0.95, na.rm = TRUE)
), by = .(Year, HTarget, Hx, Trait)]

traj_sum[, HTarget := round(HTarget, 3)]
traj_sum[, HT := factor(
  c("Low H", "Typical H", "High H")[match(HTarget, c(0.068, 0.127, 0.194))],
  levels = c("Low H", "Typical H", "High H")
)]

# ============================================================
# Compute log change, percent change, and Haldanes
# ============================================================

# Adjust generation time as appropriate for perch
gen_time_years <- 4
n_generations <- (2025 - 1900) / gen_time_years

# Extract baseline trajectories at t0 and t1
# Collapse the two stochastic replicates within each posterior particle and generation
rate_dt_u <- rate_dt[, .(
  LT = mean(LT, na.rm = TRUE),
  g = mean(g, na.rm = TRUE),
  hmax = mean(hmax, na.rm = TRUE),
  Vp_LT = mean(Vp_LT, na.rm = TRUE),
  Vp_g = mean(Vp_g, na.rm = TRUE),
  Vp_hmax = mean(Vp_hmax, na.rm = TRUE)
), by = .(Ve, sd, K, HTarget, Hx, thetaIdx, G)]

# Means wide
means_wide <- dcast(
  rate_dt_u,
  Ve + sd + K + HTarget + Hx + thetaIdx ~ G,
  value.var = c("LT", "g", "hmax")
)

# Phenotypic variances wide
vp_wide <- dcast(
  rate_dt_u,
  Ve + sd + K + HTarget + Hx + thetaIdx ~ G,
  value.var = c("Vp_LT", "Vp_g", "Vp_hmax")
)


rates <- merge(
  means_wide,
  vp_wide,
  by = c("Ve", "sd", "K", "HTarget", "Hx", "thetaIdx")
)

gen_time_years <- 4
n_generations <- (2025 - 1900) / gen_time_years

rates[, `:=`(
  sdP_LT   = sqrt((Vp_LT_10000   + Vp_LT_10125)   / 2),
  sdP_g    = sqrt((Vp_g_10000    + Vp_g_10125)    / 2),
  sdP_hmax = sqrt((Vp_hmax_10000 + Vp_hmax_10125) / 2)
)]

rates[, `:=`(
  log_change_LT   = log(LT_10125   / LT_10000),
  log_change_g    = log(g_10125    / g_10000),
  log_change_hmax = log(hmax_10125 / hmax_10000),
  
  pct_change_LT   = 100 * (LT_10125   - LT_10000)   / LT_10000,
  pct_change_g    = 100 * (g_10125    - g_10000)    / g_10000,
  pct_change_hmax = 100 * (hmax_10125 - hmax_10000) / hmax_10000,
  
  haldane_LT   = (LT_10125   - LT_10000)   / (sdP_LT   * n_generations),
  haldane_g    = (g_10125    - g_10000)    / (sdP_g    * n_generations),
  haldane_hmax = (hmax_10125 - hmax_10000) / (sdP_hmax * n_generations)
)]

# Long table for summary
rates_long <- rbindlist(list(
  rates[, .(HTarget, Hx, thetaIdx, trait = "LT",
            log_change = log_change_LT,
            pct_change = pct_change_LT,
            haldane = haldane_LT)],
  rates[, .(HTarget, Hx, thetaIdx, trait = "g",
            log_change = log_change_g,
            pct_change = pct_change_g,
            haldane = haldane_g)],
  rates[, .(HTarget, Hx, thetaIdx, trait = "hmax",
            log_change = log_change_hmax,
            pct_change = pct_change_hmax,
            haldane = haldane_hmax)]
))

haldane_summary <- rates_long[, .(
  mean_haldane = mean(haldane, na.rm = TRUE),
  q05_haldane  = quantile(haldane, 0.05, na.rm = TRUE),
  q50_haldane  = quantile(haldane, 0.50, na.rm = TRUE),
  q95_haldane  = quantile(haldane, 0.95, na.rm = TRUE),
  
  mean_pct = mean(pct_change, na.rm = TRUE),
  q05_pct  = quantile(pct_change, 0.05, na.rm = TRUE),
  q50_pct  = quantile(pct_change, 0.50, na.rm = TRUE),
  q95_pct  = quantile(pct_change, 0.95, na.rm = TRUE)
), by = .(HTarget, Hx, trait)]

haldane_summary[,plot(q50_haldane,q50_pct)]
haldane_summary[,lm(q50_pct~q50_haldane)$coeff]
