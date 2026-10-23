# ============================================================
# Figure S3: ABC--SMC diagnostics
# ============================================================
# This script produces diagnostic plots for ABC--SMC inference:
# (a) relative reduction in tolerance (epsilon),
# (b) effective sample size per accepted particle,
# (c) acceptance rate curves across generations.
# ============================================================

# ---- Packages ----------------------------------------------
library(data.table)
library(ggplot2)
library(ggh4x)        # geom_pointpath
library(viridis)     # color scales
library(egg)         # ggarrange

# ---- Load ABC results --------------------------------------
base_folder <- "/private/tmp/claude-539526166/-Users-petrikem-gitlab-REEIM/f7799891-e671-4825-bad9-a218acb19f48/scratchpad/final/res_perch6/"
APOT        <- 10125   
REEIM_res <- readRDS(file.path(base_folder, "REEIM_res.rds"))
# ---- Diagnostic plotting function --------------------------
data <- copy(REEIM_res)
plot_abc_diagnostics <- function(data, title_suffix = NULL) {
  
  # Number of accepted particles per generation
  n_accept <- sapply(data$particles, function(x) sum(x$keep))
  
  data$ecdfs <- lapply(data$particles, function(prtcl) ecdf(prtcl[,dist]))
  
  
  # --- (a) Relative improvement in epsilon ------------------
  eps <- data$epsilons
  rel_gain <- (eps[-length(eps)] - eps[-1]) / eps[-length(eps)]
  
  p_eps <- ggplot(
    data.table(Generation = 2:length(eps), Gain = rel_gain * 100),
    aes(Generation, Gain)
  ) +
    geom_pointpath() +
    geom_hline(yintercept = c(5, 1), linetype = c(2, 3)) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    labs(
      title = "a) Relative tolerance reduction",
      y = "Gain (%)",
      x = "ABC generation"
    )
  
  # --- (b) Effective sample size -----------------------------
  ESS_df <- data.table(
    Generation = seq_along(data$ESS),
    ESS_ratio  = data$ESS / n_accept
  )
  
  p_ESS <- ggplot(ESS_df, aes(Generation, ESS_ratio)) +
    geom_pointpath() +
    theme_bw() +
    theme(aspect.ratio = 1) +
    labs(
      title = "b) Effective sample size",
      y = expression(ESS / N[acc]),
      x = "ABC generation"
    )
  
  # --- (c) Acceptance rate curves ----------------------------
  eps_grid <- seq(min(eps) * 0.95, max(eps), length.out = 500)
  
  AR_df <- rbindlist(lapply(seq_along(data$ecdfs), function(g) {
    data.table(
      Generation = factor(g),
      epsilon    = eps_grid,
      AR         = data$ecdfs[[g]](eps_grid)
    )
  }))
  
  p_AR <- ggplot(
    AR_df,
    aes(epsilon, AR, color = Generation, group = Generation)
  ) +
    geom_line() +
    theme_bw() +
    theme(aspect.ratio = 1) +
    scale_color_viridis_d(option = "turbo") +
    labs(
      title = "c) Acceptance rate curves",
      y = "Acceptance probability",
      x = expression(epsilon)
    )
  
  # --- Combine panels ---------------------------------------
  ggarrange(p_eps, p_ESS, p_AR, ncol = 3)
}

# ---- Generate and save figure -------------------------------
fig_S3 <- plot_abc_diagnostics(
  REEIM_res,
  title_suffix = "ABC--SMC (5 generations)"
)

ggsave(
  filename = "MS/revision/assets/perch/ABC_diagnostics.pdf",
  plot     = fig_S3,
  width    = 17,
  height   = 5
)
