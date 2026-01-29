# ============================================================
# Table S2: Posterior parameter estimates (ABC--SMC, SLiM)
# ============================================================

# ---- Packages (minimal, required only) ---------------------
library(data.table)
library(ks)
library(kableExtra)
# ---- Load ABC results --------------------------------------
base_folder <- "./REEIM/out_REEIM/"
APOT        <- 10125   
REEIM_res <- readRDS(file.path(base_folder, "REEIM_res.rds"))
particles   <- REEIM_res$particles
last_gen    <- length(particles)

# Retain accepted particles from final generation
post <- particles[[last_gen]][keep == TRUE]

# ---- Helper functions --------------------------------------

# Weighted credible interval
weighted_ci <- function(x, w, prob = c(0.05, 0.95)) {
  ord <- order(x)
  x   <- x[ord]
  w   <- w[ord] / sum(w)
  cw  <- cumsum(w)
  approx(cw, x, prob)$y
}

# Posterior summaries
posterior_summary <- function(data, params, prob = c(0.05, 0.95)) {
  rbindlist(lapply(params, function(p) {
    values  <- data[[p]]
    weights <- data$weights
    
    mu <- sum(values * weights)
    
    ci <- weighted_ci(values, weights, prob)
    
    d  <- kde(x = values, w = weights * length(weights))
    
    md <- d$eval.points[which.max(d$estimate)]
    
    data.table(
      Parameter = p,
      Mean      = mu,
      Median    = median(d$eval.points),
      Mode      = md,
      CI_lower  = ci[1],
      CI_upper  = ci[2]
    )
  }))
}

# ---- Parameters inferred via ABC--SMC ----------------------
params <- c("g", "h", "f", "w", "Shmax", "SLT", "LT", "H")

tab <- posterior_summary(post, params,prob = c(0.05, 0.95))
#tab <- posterior_summary(post, params,prob = c(0.1, 0.9))

# ---- Format for LaTeX --------------------------------------
tab[, Mean := sprintf("%.3f [%.3f, %.3f]", Mean, CI_lower, CI_upper)]
tab[, Mode := sprintf("%.3f", Mode)]

tab[, Parameter := c(
  "$g$", "$h$", "$f$", "$w$",
  "$S_{h_{\\max}}$", "$S_{L_T}$",
  "$L_T$", "$H$"
)]

# ---- Build LaTeX table -------------------------------------
caption_text <- paste(
  "Posterior summaries of parameter estimates inferred using",
  "Approximate Bayesian Computation with Sequential Monte Carlo",
  "(ABC--SMC) and an eco-genetic individual-based model implemented",
  "in \\texttt{SLiM}. Credible intervals correspond to the central",
  "90\\% posterior mass."
)

kbl <- kable(
  tab[, .(Parameter, Mean, Mode)],
  format      = "latex",
  booktabs    = TRUE,
  escape      = FALSE,
  col.names   = c("Parameter", "Mean [95\\% CI]", "Mode"),
  caption     = caption_text,
  label       = "tab:S2_ABC_parameters"
) |>
  kable_styling(latex_options = "hold_position")

# ---- Save --------------------------------------------------
save_kable(kbl, "./Tables/Table_S2_ABC_parameters.tex")
