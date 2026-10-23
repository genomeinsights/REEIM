############################################################
##  Shared helpers for the in-silico recovery test.
##  Settings mirror REEIM/run_REEIM_ABC.R; copied rather than
##  sourced so this folder does not depend on the driver script.
############################################################

APOT      <- 10125
t0        <- 10000          # start of the fishing century (Tdelta)
BW        <- 2
from      <- 6
to        <- 40
par_names <- c("H", "f", "w", "Shmax", "SLT")
max_age   <- 26             # oldest age class in the empirical age data

get_param <- function(key, bname) {
  as.numeric(strsplit(strsplit(bname, key)[[1]][2], "_", fixed = TRUE)[[1]][1])
}

extract_one <- function(file) {
  trait_file <- sub("demographicData_", "traitStats_", file)
  if (!file.exists(trait_file)) return(NULL)
  tr <- fread(trait_file, showProgress = FALSE)
  if (max(tr$G) < APOT || !(t0 %in% tr$G)) return(NULL)

  pars <- sapply(par_names, get_param, basename(file))
  sim  <- readLines(file)
  row_at <- function(key) {
    x <- strsplit(sim[startsWith(sim, paste(key, APOT, ""))][1], " ", fixed = TRUE)[[1]]
    suppressWarnings(as.numeric(x[-(1:2)]))
  }

  ## size density at APOT (size classes are 1-cm bins, 0 = 0-1 cm)
  counts <- row_at("countSize")
  counts[is.na(counts)] <- 0
  sizes  <- rep(seq_along(counts) - 1, counts)
  dens   <- density(sizes, bw = BW, from = from, to = to)$y

  ## mean size (mm) at age 0..max_age at APOT, NA where the age class is empty
  ms     <- row_at("meanSize")
  growth <- rep(NA_real_, max_age + 1)
  n      <- min(length(ms), max_age + 1)
  growth[seq_len(n)] <- ms[seq_len(n)]
  growth[!is.finite(growth)] <- NA

  ## century trait responses (log change, as in R/trati_cor_with_H.R)
  v  <- function(col, at) tr[[col]][match(at, tr$G)]
  r  <- function(col) log(v(col, APOT)) - log(v(col, t0))
  list(
    pars   = c(pars,
               LT_APOT  = v("LT", APOT),
               r_LT     = r("LT"),
               r_g      = r("g"),
               r_hmax   = r("hmax"),
               r_h      = r("h")),
    dens   = dens,
    growth = growth,
    file   = gsub("demographicData_", "", basename(file))
  )
}


## Stack a list of extract_one() results
bind_draws <- function(res) {
  bad <- vapply(res, inherits, TRUE, "try-error")
  if (any(bad)) warning(sum(bad), " draws failed to read: ", res[bad][[1]])
  res <- Filter(function(x) is.list(x) && !inherits(x, "try-error"), res)
  list(
    theta  = as.data.table(do.call(rbind, lapply(res, `[[`, "pars"))),
    dens   = do.call(rbind, lapply(res, `[[`, "dens")),
    growth = do.call(rbind, lapply(res, `[[`, "growth")),
    file   = vapply(res, `[[`, "", "file")
  )
}

## ABC distance of simulated draws `sims` to one (pseudo-)observed
## draw `obs` (both in the bind_draws() format), as in the driver:
## SRMSE_dens + SRMSE_growth + SAD_LT. age_w is the empirical age
## composition and vw the observed within-age size variance per age
## (gc, from growth_constants()). The driver scores growth
## on individual fish, so its squared error at age a is vw[a] plus the
## squared difference in mean size; adding vw keeps the growth term on
## the same footing for pseudo-observed data, which only have means.
dist_to_target <- function(obs, sims, gc) {
  age_w <- gc$age_w; vw <- gc$vw
  d_obs <- obs$dens[1, ]
  dens  <- sqrt(rowMeans(sweep(sims$dens, 2, d_obs)^2)) / sd(d_obs)
  diff2 <- sweep(sims$growth, 2, obs$growth[1, ])^2
  w     <- matrix(age_w, nrow(diff2), ncol(diff2), byrow = TRUE) * !is.na(diff2)
  diff2 <- sweep(diff2, 2, vw, "+")
  diff2[is.na(diff2)] <- 0
  ## the driver divides by the sd of the observed sizes over the ages
  ## the simulation has, so rebuild that sd from per-age moments
  n   <- rowSums(w)
  mu  <- (w %*% gc$means0) / n
  ss  <- w %*% (vw + gc$means0^2) - n * mu^2
  sd_obs <- sqrt(ss / (n - 1))
  growth <- sqrt(rowSums(w * diff2) / n) / sd_obs
  LT    <- abs(sims$theta$LT_APOT - obs$theta$LT_APOT) / obs$theta$LT_APOT
  d <- dens + as.vector(growth) + LT
  d[!is.finite(d)] <- Inf
  d
}

## Prior support and density (uniform priors), as in the driver
out_of_bounds <- function(theta, prior_lim) {
  any(mapply(function(x, lim) x < lim[1] | x > lim[2], theta[names(prior_lim)], prior_lim))
}
prior_density <- function(theta, prior_lim) {
  prod(mapply(function(x, lim) dunif(x, lim[1], lim[2]), theta[names(prior_lim)], prior_lim))
}

## Empirical growth-scoring constants from the observed age data
growth_constants <- function(age_structure_obs, n_age) {
  a <- age_structure_obs[age < n_age]
  out <- list(age_w  = tabulate(a$age + 1, nbins = n_age),
       sd_obs = sd(age_structure_obs$size),
       vw     = vapply(0:(n_age - 1), function(x) {
         v <- a$size[a$age == x]; if (length(v)) mean((v - mean(v))^2) else 0 }, 0),
       means  = vapply(0:(n_age - 1), function(x) {
         v <- a$size[a$age == x]; if (length(v)) mean(v) else NA_real_ }, 0))
  out$means0 <- ifelse(is.na(out$means), 0, out$means)
  out
}
