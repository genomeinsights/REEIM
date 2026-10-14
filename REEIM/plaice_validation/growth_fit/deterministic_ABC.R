## =============================================================================
## Deterministic growth-model ABC-SMC (REEIM), any study system.
## Copy of deterministic_ABC/deterministic_ABC.R that can read a system's input folder
## (issue I-19). With no input folder it uses the perch settings and reproduces the
## original script.
##
## Fits deterministic growth curves to individual length-at-age data with
## ABC-SMC and compares five models:
##   (1) von Bertalanffy                         (also fitted by nls)
##   (2) Lester et al. (2004): h0 = 1, g_i = 0, L_T fixed
##   (3) full tri-phasic model: h, h0, g_i, L_T, g all free
##   (4) reduced model: g_i = 0, L_T = 0 fixed (one g for all ages after year 1)
##   (5) tri-phasic model with L_T fixed (the model in the Methods)
##   (6) as (5) with g also fixed (only if growth_priors.csv gives g_fixed; issue I-18)
##
## Usage:
##   Rscript deterministic_ABC.R [observed_data.RData] [cores] [input folder] [output folder]
## Defaults: ../empirical_data/observed_data.RData, all cores but one, no input folder
## (perch settings), ./output. With an input folder:
##   L_T     = LT_obs in <input>/abc_settings.csv
##   max_age = maxAge + 1 in <input>/parameters.csv (perch: 25 + 1 = 26, as before)
##   priors  = <input>/growth_priors.csv (parameter,lower,upper), overriding the perch
##             bounds below for the parameters it lists; a row g_fixed (lower = upper)
##             adds model (6).
## Example (plaice, from the repository root):
##   Rscript REEIM/plaice_validation/growth_fit/deterministic_ABC.R \
##     REEIM/plaice_validation/data/observed_data_plaice_female.RData 4 \
##     REEIM/input/plaice REEIM/plaice_validation/growth_fit/output
## The .RData file must contain a data.table `age_structure_obs` with columns
## `age` (years) and `size` (mm), one row per aged fish.
## Outputs are written to the output folder.
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ks)        # weighted KDE for posterior modes
  library(mvtnorm)   # Gaussian-mixture proposal density
  library(ggplot2)
  library(parallel)
})

args      <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[1] else "../empirical_data/observed_data.RData"
cores     <- if (length(args) >= 2) as.integer(args[2]) else max(1L, detectCores() - 1L)
input_dir <- if (length(args) >= 3) args[3] else NA_character_
out_dir   <- if (length(args) >= 4) args[4] else "output"
if (.Platform$OS.type == "windows") cores <- 1L   # mclapply forks only on Unix

## ---- Settings (as in the original analysis) --------------------------------
N        <- 2500                 # particles simulated per generation
Gens     <- 5                    # SMC generations
alpha    <- 0.2                  # fraction of particles accepted per generation
c_cov    <- 2                    # covariance inflation of the proposal kernel
ci_prob  <- c(0.05, 0.95)        # 90% credible interval, as in Table S2
max_age  <- 26
L_T_ref  <- 169.2                # literature size at maturation (mm)

## Perch prior bounds (the published analysis); growth_priors.csv can override them.
bounds <- list(linf = c(400, 500), k = c(0.02, 0.12), t0 = c(-4, 0),
               h = c(10, 100), h0 = c(1, 3), g_i = c(0, 0.5), g = c(0, 0.5), L_T = c(0, 200))
g_fixed <- NA_real_

read_settings <- function(file) {
  d <- read.csv(file, colClasses = "character", comment.char = "#", strip.white = TRUE)
  setNames(as.list(d$value), d$name)
}
if (!is.na(input_dir)) {
  L_T_ref <- as.numeric(read_settings(file.path(input_dir, "abc_settings.csv"))$LT_obs)
  max_age <- as.integer(read_settings(file.path(input_dir, "parameters.csv"))$maxAge) + 1L
  pf <- file.path(input_dir, "growth_priors.csv")
  if (file.exists(pf)) {
    gp <- read.csv(pf, comment.char = "#", strip.white = TRUE)
    for (i in seq_len(nrow(gp))) {
      if (gp$parameter[i] == "g_fixed") g_fixed <- gp$lower[i]
      else bounds[[gp$parameter[i]]] <- c(gp$lower[i], gp$upper[i])
    }
  }
  cat("Input folder", input_dir, ": L_T =", L_T_ref, " max_age =", max_age, " g_fixed =", g_fixed, "\n")
}
B <- function(...) bounds[c(...)]
set.seed(20251014)
RNGkind("L'Ecuyer-CMRG")         # reproducible streams under mclapply

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
load(data_file)
stopifnot(exists("age_structure_obs"), all(c("age", "size") %in% names(age_structure_obs)))
obs <- as.data.table(age_structure_obs)[, .(age, size)]

## ---- Growth models -----------------------------------------------------------
## Each returns predicted length (mm) at ages 0:max_age for a named vector theta.

growth_vB <- function(theta, max_age) {
  age <- 0:max_age
  theta[["linf"]] * (1 - exp(-theta[["k"]] * (age - theta[["t0"]])))
}

## Tri-phasic model (manuscript Methods, "Growth"):
##   first year:            L(0) = h * h0
##   immature (L <  L_T):   L(t+1) = 3/(g_i + 3) * (L(t) + h)
##   mature   (L >= L_T):   L(t+1) = 3/(g   + 3) * (L(t) + h)
## which is dL = h - g/(g+3) * (L + h). Maturity is absorbing.
growth_triphasic <- function(theta, max_age) {
  L <- numeric(max_age + 1)
  L[1] <- theta[["h"]] * theta[["h0"]]
  mature <- FALSE
  for (t in 1:max_age) {
    mature <- mature || L[t] >= theta[["L_T"]]
    g_t <- if (mature) theta[["g"]] else theta[["g_i"]]
    L[t + 1] <- 3 / (g_t + 3) * (L[t] + theta[["h"]])
  }
  L
}

## Standardised RMSE between predicted and observed individual lengths.
distance <- function(theta, curve_fun, obs, max_age) {
  pred <- curve_fun(theta, max_age)[obs$age + 1]
  sqrt(mean((pred - obs$size)^2)) / sd(obs$size)
}

## ---- Priors --------------------------------------------------------------------
## `free`: uniform prior bounds; `fixed`: parameters held constant.
models <- list(
  vB = list(
    label = "(1) von Bertalanffy", curve = growth_vB,
    free  = B("linf", "k", "t0"),
    fixed = list()),
  Lester04 = list(
    label = "(2) Lester et al. 2004", curve = growth_triphasic,
    free  = B("h", "g"),
    fixed = list(h0 = 1, g_i = 0, L_T = L_T_ref)),
  full = list(
    label = "(3) Full tri-phasic", curve = growth_triphasic,
    free  = B("h", "h0", "g_i", "L_T", "g"),
    fixed = list()),
  reduced = list(
    label = "(4) Reduced (single g)", curve = growth_triphasic,
    free  = B("h", "h0", "g"),
    fixed = list(g_i = 0, L_T = 0)),
  triphasic_LT = list(
    label = "(5) Tri-phasic, L_T fixed", curve = growth_triphasic,
    free  = B("h", "h0", "g_i", "g"),
    fixed = list(L_T = L_T_ref))
)
if (!is.na(g_fixed))
  models$triphasic_LT_g <- list(
    label = sprintf("(6) Tri-phasic, L_T and g = %g fixed", g_fixed), curve = growth_triphasic,
    free  = B("h", "h0", "g_i"),
    fixed = list(L_T = L_T_ref, g = g_fixed))

## ---- ABC-SMC -----------------------------------------------------------------------
## Generation 1 samples the prior. Later generations draw from a weighted
## Gaussian mixture centred on the previous generation's accepted particles with
## covariance c_cov * (weighted covariance); proposals outside the prior are
## redrawn. The alpha-quantile of each generation's distances sets epsilon, and
## importance weights are w_i ∝ p(theta_i) / q_t(theta_i), where q_t is that
## same mixture, so sampling and weighting use one proposal density.

in_bounds <- function(x, lo, hi) all(x >= lo & x <= hi)

abc_smc <- function(model, obs, N, Gens, alpha, c_cov, max_age, cores) {
  lo <- sapply(model$free, `[`, 1)
  hi <- sapply(model$free, `[`, 2)
  par_names <- names(model$free)
  full_theta <- function(x) c(setNames(x, par_names), unlist(model$fixed))
  dist_fun <- function(x) distance(full_theta(x), model$curve, obs, max_age)
  log_prior <- -sum(log(hi - lo))      # uniform prior: constant inside bounds

  generations <- vector("list", Gens)
  diagnostics <- data.table(generation = 1:Gens, epsilon = NA_real_, ESS = NA_real_,
                            in_bounds_rate = NA_real_)
  prev <- NULL
  for (GEN in 1:Gens) {
    cat(model$label, "- generation", GEN, "\n")
    if (GEN == 1) {
      draw <- function(i) {
        x <- runif(length(lo), lo, hi)
        c(x, dist = dist_fun(x), tries = 1)
      }
    } else {
      Sigma <- c_cov * cov.wt(prev$theta, wt = prev$w)$cov
      draw <- function(i) {
        tries <- 0
        repeat {
          tries <- tries + 1
          j <- sample.int(nrow(prev$theta), 1, prob = prev$w)
          x <- as.vector(rmvnorm(1, prev$theta[j, ], Sigma))
          if (in_bounds(x, lo, hi)) break
        }
        c(x, dist = dist_fun(x), tries = tries)
      }
    }
    sims <- do.call(rbind, mclapply(1:N, draw, mc.cores = cores))
    colnames(sims) <- c(par_names, "dist", "tries")

    eps <- quantile(sims[, "dist"], alpha)
    acc <- sims[sims[, "dist"] < eps, , drop = FALSE]
    theta <- acc[, par_names, drop = FALSE]

    if (GEN == 1) {
      w <- rep(1, nrow(theta))
    } else {
      w <- unlist(mclapply(1:nrow(theta), function(i) {
        q <- sum(prev$w * dmvnorm(prev$theta, mean = theta[i, ], sigma = Sigma))
        exp(log_prior) / q
      }, mc.cores = cores))
    }
    w <- w / sum(w)

    diagnostics[GEN, `:=`(epsilon = eps, ESS = 1 / sum(w^2),
                          in_bounds_rate = N / sum(sims[, "tries"]))]
    generations[[GEN]] <- cbind(theta, dist = acc[, "dist"], weights = w)
    prev <- list(theta = theta, w = w)
  }
  names(generations) <- paste("Generation", 1:Gens)
  list(particles = generations, diagnostics = diagnostics, model = model)
}

## ---- Posterior summaries ------------------------------------------------------------
weighted_ci <- function(samples, weights, prob) {
  o <- order(samples)
  approx(cumsum(weights[o]), samples[o], xout = prob, rule = 2, ties = "ordered")$y
}

posterior_summary <- function(res, prob = ci_prob) {
  final <- res$particles[[length(res$particles)]]
  w <- final[, "weights"]
  free <- rbindlist(lapply(names(res$model$free), function(p) {
    x <- final[, p]
    d <- kde(x = x, w = w * length(w))
    ci <- weighted_ci(x, w, prob)
    data.table(Parameter = p, Mean = sum(x * w), Mode = d$eval.points[which.max(d$estimate)],
               CI_lower = ci[1], CI_upper = ci[2],
               Prior = sprintf("U(%g, %g)", res$model$free[[p]][1], res$model$free[[p]][2]))
  }))
  fixed <- if (length(res$model$fixed)) data.table(
    Parameter = names(res$model$fixed), Mean = unlist(res$model$fixed),
    Mode = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_, Prior = "fixed") else NULL
  cbind(Model = res$model$label, rbind(free, fixed))
}

## Warn when a credible interval reaches a prior bound (prior may be truncating).
check_prior_bounds <- function(summ, model, tol = 0.01) {
  for (p in names(model$free)) {
    b <- model$free[[p]]; s <- summ[Parameter == p]
    span <- diff(b)
    if (s$CI_lower - b[1] < tol * span || b[2] - s$CI_upper < tol * span)
      warning(sprintf("%s: credible interval for %s touches the prior bound U(%g, %g)",
                      model$label, p, b[1], b[2]), call. = FALSE)
  }
}

## ---- Run ------------------------------------------------------------------------------
results <- lapply(models, abc_smc, obs = obs, N = N, Gens = Gens, alpha = alpha,
                  c_cov = c_cov, max_age = max_age, cores = cores)
saveRDS(results, file.path(out_dir, "ABC_results.rds"))

summaries <- rbindlist(lapply(results, posterior_summary))
for (m in names(models)) check_prior_bounds(summaries[Model == models[[m]]$label], models[[m]])

## von Bertalanffy by nonlinear least squares, for reference
vB_nls <- nls(size ~ linf * (1 - exp(-k * (age - t0))), data = obs,
              start = list(linf = 450, k = 0.08, t0 = -1))
vB_nls_coef <- coef(vB_nls)
vB_nls_ci <- confint.default(vB_nls, level = diff(ci_prob))
summaries <- rbind(summaries, data.table(
  Model = "(1) von Bertalanffy, nls", Parameter = names(vB_nls_coef), Mean = vB_nls_coef,
  Mode = NA_real_, CI_lower = vB_nls_ci[, 1], CI_upper = vB_nls_ci[, 2], Prior = "-"))

diagnostics <- rbindlist(lapply(results, function(r) cbind(Model = r$model$label, r$diagnostics)))
fwrite(summaries, file.path(out_dir, "posterior_summaries.csv"))
fwrite(diagnostics, file.path(out_dir, "ABC_diagnostics.csv"))
print(summaries[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])
print(diagnostics)

if (requireNamespace("knitr", quietly = TRUE)) {
  fmt <- summaries[, lapply(.SD, function(x) if (is.numeric(x)) formatC(x, digits = 3, format = "fg") else x)]
  writeLines(knitr::kable(fmt, format = "latex", booktabs = TRUE,
                          caption = "ABC-SMC posterior summaries of deterministic growth models"),
             file.path(out_dir, "posterior_summaries.tex"))
}

## ---- Figures ------------------------------------------------------------------------------
model_levels <- c("(1) von Bertalanffy, nls", sapply(models, `[[`, "label"))
full_theta <- function(res, x) c(x[names(res$model$free)], unlist(res$model$fixed))

## Curves at the posterior mean of each model
curves <- rbindlist(c(
  list(data.table(age = 0:max_age, size = growth_vB(vB_nls_coef, max_age),
                  Model = "(1) von Bertalanffy, nls")),
  lapply(results, function(r) {
    s <- posterior_summary(r)
    theta <- setNames(s$Mean, s$Parameter)
    data.table(age = 0:max_age, size = r$model$curve(theta, max_age), Model = r$model$label)
  })))
curves[, Model := factor(Model, levels = model_levels)]

p_curves <- ggplot(obs, aes(age, size)) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 0.8) +
  geom_line(data = curves, aes(col = Model), linewidth = 0.8) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Age (years)", y = "Length (mm)", title = "Growth curves fitted to perch data") +
  theme_bw() + theme(aspect.ratio = 0.5)
ggsave(file.path(out_dir, "growth_curves_posterior_mean.pdf"), p_curves, width = 8, height = 4.5)

## Final-generation particle curves for each ABC model
particle_curves <- rbindlist(lapply(results, function(r) {
  final <- r$particles[[Gens]]
  rbindlist(lapply(seq_len(nrow(final)), function(i) {
    data.table(particle = i, age = 0:max_age,
               size = r$model$curve(full_theta(r, final[i, ]), max_age), Model = r$model$label)
  }))
}))
particle_curves[, Model := factor(Model, levels = model_levels)]
p_particles <- ggplot(obs, aes(age, size, group = age)) +
  geom_boxplot(alpha = 0.2, outlier.size = 0.5) +
  geom_line(data = particle_curves, aes(age, size, group = particle),
            col = "steelblue", alpha = 0.05, inherit.aes = FALSE) +
  geom_line(data = curves[Model != "(1) von Bertalanffy, nls"], aes(age, size),
            col = "black", linewidth = 1, inherit.aes = FALSE) +
  facet_wrap(~Model, ncol = 2) +
  labs(x = "Age (years)", y = "Length (mm)") +
  theme_bw()
ggsave(file.path(out_dir, "growth_curves_particles.pdf"), p_particles, width = 9, height = 9)

## Distances of accepted particles in the final generation
dists <- rbindlist(lapply(results, function(r)
  data.table(Model = r$model$label, dist = r$particles[[Gens]][, "dist"])))
dists[, Model := factor(Model, levels = model_levels)]
p_dist <- ggplot(dists, aes(dist, col = Model, fill = Model)) +
  geom_density(alpha = 0.5) +
  scale_color_brewer(palette = "Dark2", drop = FALSE) +
  scale_fill_brewer(palette = "Dark2", drop = FALSE) +
  labs(x = expression(italic(rho(theta))), y = "Density",
       title = "Distances between observed and simulated data") +
  theme_bw() + theme(aspect.ratio = 0.5)
ggsave(file.path(out_dir, "distance_distributions.pdf"), p_dist, width = 8, height = 4.5)

cat("Done. Results in", normalizePath(out_dir), "\n")
