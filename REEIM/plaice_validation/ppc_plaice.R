## Posterior predictive check (PPC) of evolutionary change for a plaice ABC: can the fitted model produce the
## published female age-4 Lp50 (PMRN) declines?
##
## Usage (from the repository root):
##   Rscript REEIM/plaice_validation/ppc_plaice.R genetic <abc_out> <ppc_out>
##       Weighted change in genetic mean L_T over the accepted particles of the last ABC generation
##       (uses their traitStats files; no new simulations). Also a naive evolutionary sensitivity, change / (10 H)
##       in % per 0.1 of H, which assumes no change without fishing (see the paired H = 0 controls below).
##   Rscript REEIM/plaice_validation/ppc_plaice.R prepare <abc_out> <ppc_out> <input_dir> [n=250] [seed=1]
##       Draws n parameter sets from the weighted posterior and writes <ppc_out>/design.csv and cmds.txt for
##       REEIM_vesurvey.slim (REEIM.slim plus a random survey of 2000 females every year, 1900-2020). Run cmds.txt with
##       GNU parallel, e.g.  cat <ppc_out>/cmds.txt | parallel -j <cores>
##       Each draw i also gets a paired control run (rep n + i) with H = 0 and the same seed: the burn-in is
##       identical, so the difference between the pair is the effect of fishing on L_T (drift diverges after 1850).
##   Rscript REEIM/plaice_validation/ppc_plaice.R survey <ppc_out>
##       Fits ogive L50 (per survey year) and PMRN Lp50 at age 4 (per decade) to each simulated survey, and summarises
##       genetic L_T, ogive and PMRN changes. The PMRN is estimated two ways: with first-time maturing known (logistic
##       first ~ age + length; not observable in real surveys) and with the Barot et al. (2004) demographic method
##       from age-specific ogives and growth, as in the published plaice PMRNs (Grift et al. 2003; van Walraven et al. 2010).
##       Evolutionary sensitivity: per draw, (change in genetic L_T at the drawn H - change in its H = 0 control) / (10 H),
##       i.e. the fishing-induced change in L_T in % per 0.1 of terminal harvest rate, per window.
## Published (female, age 4): Grift et al. 2003 Lp50 cohorts 1955-64 -> 1985-94: 34.4 -> 30.1 cm (-12.5%);
## van Walraven et al. 2010 Lp50 1900s -> 1980s -> 2000s: 41.6 -> 33.3 -> 22.9 cm (-20%, then -31%).
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
args <- commandArgs(TRUE); mode <- args[1]
Td <- 10000; start_year <- 1850
win <- list(`1905-85` = c(1905, 1985), `1960-90` = c(1960, 1990), `1985-2005` = c(1985, 2005), `1900-2020` = c(1900, 2020))
pub <- c(`1905-85` = -20.0, `1960-90` = -12.5, `1985-2005` = -31.2, `1900-2020` = NA)
wq <- function(x, w, p) { o <- order(x); x[o][which(cumsum(w[o]) / sum(w) >= p)[1]] }
summ <- function(d, w = rep(1, nrow(d)), lab) rbindlist(lapply(names(win), function(k) {
  x <- d[[k]]; ok <- !is.na(x); x <- x[ok]; ww <- w[ok]
  data.table(measure = lab, window = k, n = length(x), median = wq(x, ww, 0.5), q05 = wq(x, ww, 0.05), q95 = wq(x, ww, 0.95),
             published = pub[k], P_reach_published = if (is.na(pub[k])) NA_real_ else sum(ww[x <= pub[k]]) / sum(ww))
}))
change <- function(tr, y) sapply(win, function(v) 100 * (tr[year == v[2], value][1] / tr[year == v[1], value][1] - 1))

if (mode == "genetic") {
  abc_out <- args[2]; out <- args[3]; dir.create(out, showWarnings = FALSE, recursive = TRUE)
  r <- readRDS(file.path(abc_out, "REEIM_res.rds")); G <- length(r$particles); p <- r$particles[[G]][keep == TRUE]
  ch <- rbindlist(lapply(seq_len(nrow(p)), function(i) {
    f <- file.path(abc_out, paste0("GEN_", G), paste0("traitStats_", p$file[i])); if (!file.exists(f)) return(NULL)
    tr <- fread(f)[, .(year = G - Td + start_year, value = LT)]
    data.table(i = i, weight = p$weights[i], H = p$H[i], as.data.table(as.list(change(tr))))
  }))
  setnames(ch, names(ch), sub("\\.", "-", names(ch)))
  sens <- copy(ch)[, (names(win)) := lapply(.SD, function(x) x / (10 * H)), .SDcols = names(win)]
  s <- rbind(summ(ch, ch$weight, "genetic mean L_T (accepted particles)"),
             summ(sens, sens$weight, "naive sensitivity: % change in L_T per 0.1 H (accepted particles)")[, `:=`(published = NA, P_reach_published = NA)])[
    , `:=`(generation = G, epsilon = r$epsilons[G])]
  fwrite(ch, file.path(out, "ppc_genetic_particles.csv")); fwrite(s, file.path(out, "ppc_genetic_summary.csv"))
  print(s, digits = 3)
}

if (mode == "prepare") {
  abc_out <- args[2]; dir.create(args[3], showWarnings = FALSE, recursive = TRUE); out <- normalizePath(args[3]); input_dir <- normalizePath(args[4])
  n <- if (length(args) >= 5) as.integer(args[5]) else 250; set.seed(if (length(args) >= 6) as.integer(args[6]) else 1)
  slim_script <- normalizePath("REEIM/plaice_validation/test_runs/REEIM_vesurvey.slim")
  r <- readRDS(file.path(abc_out, "REEIM_res.rds")); G <- length(r$particles); p <- r$particles[[G]][keep == TRUE]
  pars <- fread(file.path(input_dir, "priors.csv"))$parameter
  idx <- sample(nrow(p), n, replace = TRUE, prob = p$weights)   # draws from the weighted posterior (no kernel jitter)
  d <- data.table(set = "ppc", rep = seq_len(n), p[idx, ..pars], particle = idx, control_for = NA_integer_)
  d0 <- copy(d)[, `:=`(set = "ppc_H0", control_for = rep, rep = rep + n, H = 0)]   # paired no-fishing controls
  d <- rbind(d, d0)
  dir.create(file.path(out, "out"), showWarnings = FALSE, recursive = TRUE)
  fwrite(d, file.path(out, "design.csv"))
  cmds <- sapply(seq_len(2 * n), function(i) paste0("slim -s ", 7000 + ((i - 1) %% n) + 1, " -d \"folder='", file.path(out, "out"), "'\" -d \"inputDir='", input_dir,
    "'\" -d rep=", d$rep[i], " ", paste0("-d ", pars, "=", as.character(signif(unlist(d[i, ..pars]), 4)), collapse = " "), " ", slim_script, " > /dev/null 2>&1"))
  writeLines(cmds, file.path(out, "cmds.txt")); cat("Wrote", 2 * n, "commands (", n, "draws +", n, "H = 0 controls) to", file.path(out, "cmds.txt"), "\n")
}

if (mode == "survey") {
  out <- args[2]; dall <- fread(file.path(out, "design.csv")); dir.create(file.path(out, "fig"), showWarnings = FALSE)
  d <- if ("control_for" %in% names(dall)) dall[is.na(control_for)] else dall
  l50 <- function(x) { m <- glm(mature ~ length, binomial, x); b <- coef(m); -b[1] / b[2] }
  lp50 <- function(x, a = 4) { if (sum(x$first) < 10) return(NA_real_); b <- coef(glm(first ~ age + length, binomial, x)); -(b[1] + b[2] * a) / b[3] }
  ## Barot et al. (2004) demographic method, as used for the published plaice PMRNs: from the maturity ogives at
  ## ages a-1 and a and the mean growth increment ds between them (first-time maturing is not observed),
  ## p(a, s) = [o(a, s) - o(a-1, s - ds)] / [1 - o(a-1, s - ds)]; Lp50 = length where p(a, s) = 0.5.
  ## Ogives: one logistic over ages a-1 and a with an age effect and a common length slope (separate per-age
  ## slopes often give p(a, s) curves that never reach 0.5).
  barot <- function(x, a = 4) {
    x <- x[age %in% c(a - 1, a)]; if (nrow(x) < 40 || length(unique(x$mature)) < 2) return(NA_real_)
    m <- suppressWarnings(glm(mature ~ factor(age) + length, binomial, x)); ds <- x[age == a, mean(length)] - x[age == a - 1, mean(length)]
    p <- function(L) { o1 <- predict(m, data.frame(age = a - 1, length = L - ds), type = "response")
      o2 <- predict(m, data.frame(age = a, length = L), type = "response"); (o2 - o1) / (1 - o1) }
    r <- try(uniroot(function(L) p(L) - 0.5, c(50, 1000))$root, silent = TRUE); if (inherits(r, "try-error")) NA_real_ else r
  }
  res <- list(); traj <- list()
  for (i in seq_len(nrow(d))) {
    sf <- Sys.glob(sprintf("%s/out/survey_*_rep%d.txt", out, d$rep[i])); tf <- Sys.glob(sprintf("%s/out/traitStats_*_rep%d.txt", out, d$rep[i]))
    if (!length(sf) || !length(tf)) next
    tr <- fread(tf)[, .(year = G - Td + start_year, value = LT)]; if (max(tr$year) < 2020) next
    s <- fread(sf); s[, decade := 10 * (year %/% 10)]
    og <- s[, .(value = l50(.SD)), by = year]
    pm <- s[mature == 0 | first == 1, .(value = lp50(.SD)), by = decade][, .(year = decade + 5, value)]
    pb <- s[, .(value = barot(.SD)), by = decade][, .(year = decade + 5, value)]
    ## PMRN windows use decade midpoints: 1905-1985, 1965-1995 (~1960-90), 1985-2005, 1905-2015
    dch <- function(pm) { v <- function(a, b) { x <- 100 * (pm[year == b, value] / pm[year == a, value] - 1); if (length(x)) x else NA_real_ }
      c(`1905-85` = v(1905, 1985), `1960-90` = v(1965, 1995), `1985-2005` = v(1985, 2005), `1900-2020` = v(1905, 2015)) }
    res[[length(res) + 1]] <- rbind(data.table(rep = d$rep[i], measure = "genetic mean L_T", as.data.table(as.list(change(tr)))),
                                    data.table(rep = d$rep[i], measure = "ogive L50 (survey)", as.data.table(as.list(change(og)))),
                                    data.table(rep = d$rep[i], measure = "PMRN Lp50 age 4 (survey, first-time maturing known)", as.data.table(as.list(dch(pm)))),
                                    data.table(rep = d$rep[i], measure = "PMRN Lp50 age 4 (survey, Barot 2004 method)", as.data.table(as.list(dch(pb)))), use.names = FALSE)
    ref <- function(x) x[, rel := 100 * (value / mean(value[year >= 1895 & year <= 1915], na.rm = TRUE) - 1)]
    traj[[length(traj) + 1]] <- rbind(ref(tr[year >= 1895 & year %% 5 == 0])[, measure := "genetic mean L_T"],
                                      ref(og)[, measure := "ogive L50 (survey)"], ref(pm)[, measure := "PMRN Lp50 age 4 (survey, first-time maturing known)"],
                                      ref(pb)[, measure := "PMRN Lp50 age 4 (survey, Barot 2004 method)"])[, rep := d$rep[i]]
  }
  res <- rbindlist(res); setnames(res, names(res), sub("\\.", "-", names(res))); traj <- rbindlist(traj)
  s <- rbindlist(lapply(unique(res$measure), function(m) summ(res[measure == m], lab = m)))
  fwrite(res, file.path(out, "ppc_survey_draws.csv")); fwrite(s, file.path(out, "ppc_survey_summary.csv")); print(s, digits = 3)
  band <- traj[, .(med = median(rel, na.rm = TRUE), lo = quantile(rel, 0.05, na.rm = TRUE), hi = quantile(rel, 0.95, na.rm = TRUE)), by = .(measure, year)]
  FI <- fread("REEIM/input/plaice/fishing_intensity.csv")   # run from the repository root
  pts <- data.table(year = c(1905, 1985, 2005, 1960, 1990), rel = c(0, -20.0, -45.0, 0, -12.5),
                    src = c(rep("van Walraven 2010", 3), rep("Grift 2003 (0 at 1955-64)", 2)))
  lo <- max(-60, min(c(band$lo, pts$rel), na.rm = TRUE)); hi <- min(30, max(c(band$hi, 5), na.rm = TRUE)); sc <- hi - lo
  FI[, y := lo + intensity * sc]
  p <- ggplot() + geom_ribbon(data = FI, aes(year, ymin = lo, ymax = y), fill = "grey90") + geom_hline(yintercept = 0, linetype = 3) +
    geom_ribbon(data = band, aes(year, ymin = lo, ymax = hi, fill = measure), alpha = 0.25) + geom_line(data = band, aes(year, med, colour = measure)) +
    geom_line(data = pts, aes(year, rel, group = src, linetype = src)) + geom_point(data = pts, aes(year, rel, shape = src)) +
    scale_y_continuous(name = "Change since 1895-1915 (%)", sec.axis = sec_axis(~ (. - lo) / sc, name = "Relative fishing intensity (grey)")) +
    coord_cartesian(xlim = c(1895, 2022), ylim = c(lo, hi), expand = FALSE) + theme_bw(base_size = 10) +
    theme(legend.position = "bottom", legend.box = "vertical") + labs(x = "Year", colour = NULL, fill = NULL, linetype = NULL, shape = NULL,
      title = sprintf("Posterior predictive: genetic L_T and survey-estimated maturation size (%d draws; median, 90%% band)", nrow(d)))
  ggsave(file.path(out, "fig/ppc_trajectories.png"), p, width = 9, height = 5.5, dpi = 150)

  ## Evolutionary sensitivity from the paired H = 0 controls
  if ("control_for" %in% names(dall) && any(!is.na(dall$control_for))) {
    gch <- function(r) { tf <- Sys.glob(sprintf("%s/out/traitStats_*_rep%d.txt", out, r)); if (!length(tf)) return(NULL)
      tr <- fread(tf)[, .(year = G - Td + start_year, value = LT)]; if (max(tr$year) < 2020) return(NULL); change(tr) }
    sv <- rbindlist(lapply(dall[!is.na(control_for), which = TRUE], function(j) {
      c0 <- gch(dall$rep[j]); m <- dall[rep == dall$control_for[j]]; c1 <- gch(m$rep)
      if (is.null(c0) || is.null(c1)) return(NULL)
      rbind(data.table(rep = m$rep, H = m$H, measure = "change at drawn H", as.data.table(as.list(c1))),
            data.table(rep = m$rep, H = m$H, measure = "change without fishing (H = 0 control)", as.data.table(as.list(c0))),
            data.table(rep = m$rep, H = m$H, measure = "fishing effect (drawn H - control)", as.data.table(as.list(c1 - c0))),
            data.table(rep = m$rep, H = m$H, measure = "sensitivity: % change in L_T per 0.1 H", as.data.table(as.list((c1 - c0) / (10 * m$H)))))
    }))
    setnames(sv, names(sv), sub("\\.", "-", names(sv)))
    ss <- rbindlist(lapply(unique(sv$measure), function(m) summ(sv[measure == m], lab = paste("genetic L_T:", m))))
    ss[!grepl("at drawn H", measure), `:=`(published = NA, P_reach_published = NA)]
    fwrite(merge(sv, d[, c("rep", setdiff(names(d), c("set", "rep", "particle", "control_for"))), with = FALSE], by = c("rep", "H")),
           file.path(out, "ppc_sensitivity_draws.csv"))
    fwrite(ss, file.path(out, "ppc_sensitivity_summary.csv")); print(ss, digits = 3)
    x <- merge(sv[measure == "sensitivity: % change in L_T per 0.1 H"], d, by = c("rep", "H"))
    pars <- intersect(c("w", "Shmax", "SLT", "f", "H"), names(x))
    cat("\nSpearman correlation of the 1900-2020 sensitivity with the parameters:\n")
    print(round(sapply(pars, function(v) cor(x[["1900-2020"]], x[[v]], method = "s", use = "complete")), 2))
  }
}
