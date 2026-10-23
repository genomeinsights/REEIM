# ============================================================
# Table S2: Posterior parameter estimates (ABC--SMC, SLiM)
# ============================================================

# ---- Packages (minimal, required only) ---------------------
library(data.table)
library(ks)
library(kableExtra)
# ---- Load ABC results --------------------------------------
base_folder <- Sys.getenv("REEIM_ABC_OUT", "/private/tmp/claude-539526166/-Users-petrikem-gitlab-REEIM/f7799891-e671-4825-bad9-a218acb19f48/scratchpad/final/res_plaice/")   # folder with REEIM_res.rds
input_dir   <- Sys.getenv("REEIM_INPUT", "./REEIM/input/perch")     # for the H prior (priors.csv)
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
      Median    = weighted_ci(values, weights, 0.5),
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

# ---- Evidence for fishing (H > 0) ---------------------------
# Savage-Dickey Bayes factor for H = 0 against H > 0 under the uniform prior:
# BF01 = posterior density at 0 / prior density at 0, estimated as the ratio of
# posterior to prior mass in [0, delta] (two deltas as a sensitivity check).
# BF01 > 1 favours no fishing; BF10 = 1 / BF01.
H_prior <- unlist(fread(file.path(input_dir, "priors.csv"))[parameter == "H", .(lower, upper)])
w_post  <- post$weights / sum(post$weights)
H_evidence <- rbindlist(lapply(c(0.005, 0.01), function(delta) {
  post_mass  <- sum(w_post[post$H < delta])
  prior_mass <- (delta - H_prior[1]) / (H_prior[2] - H_prior[1])
  data.table(delta = delta, post_mass = post_mass, prior_mass = prior_mass,
             n_particles = sum(post$H < delta),
             BF01 = post_mass / prior_mass, BF10 = prior_mass / post_mass,
             P_H_above = 1 - post_mass)
}))
## (output to the plaice figure folder)
fwrite(H_evidence, "./compare_old_new/figures_dardel_plaice/Table_S2_H_evidence_plaice.csv")
print(H_evidence)

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
  "90\\% posterior mass.",
  sprintf(paste("Evidence for fishing: the posterior probability that $H > 0.01$ is %.2f,",
                "and the Savage--Dickey Bayes factor for $H > 0$ against $H = 0$ is %.1f",
                "(posterior against prior mass below $H = 0.01$; %.1f below $H = 0.005$)."),
          H_evidence[delta == 0.01, P_H_above], H_evidence[delta == 0.01, BF10], H_evidence[delta == 0.005, BF10])
)

kbl <- kable(
  tab[, .(Parameter, Mean, Mode)],
  format      = "latex",
  booktabs    = TRUE,
  escape      = FALSE,
  col.names   = c("Parameter", "Mean [90\\% CI]", "Mode"),
  caption     = caption_text,
  label       = "tab:S2_ABC_parameters"
) |>
  kable_styling(latex_options = "hold_position")

# ---- Save --------------------------------------------------
save_kable(kbl, "./compare_old_new/figures_dardel_plaice/Table_S2_ABC_parameters_plaice.tex")
