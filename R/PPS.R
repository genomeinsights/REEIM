library(performance)
library(knitr)
library(kableExtra)
library(ggsci)
library(scales)

# Core
library(data.table)
library(lmerTest)
library(lmtest)
library(effectsize)
library(parallel)
library(ggplot2)

# Modeling
library(lme4)        # mixed-effects
library(broom.mixed) # tidy summaries (optional)
library(patchwork)
run_root <- "./REEIM/scenarios/"
files <- list.files(run_root, recursive = TRUE, full.names = TRUE, pattern = "traitStats")

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

dt_median[ Trait ==  "A  Annual growth increment", Trait2:='bold(A)~"Annual growth increment"']
dt_median[ Trait ==  "B  Size at maturation", Trait2:='bold(B)~"Size at maturation"']
dt_median[ Trait ==  "C  Parental investment", Trait2:='bold(C)~"Parental investment"']


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
  facet_wrap(~Trait2, scales = "free_y", labeller = label_parsed) +
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

#### === violin plots of scenarios === ###


run_root <- "./REEIM/scenarios/"      # folder that contains scenario subfolders
t0 <- 10000; t1 <- 10125                      # century window
trait_names <- c("LT","g","hmax")             # traits to summarize

# list all trait_stats files
assemble_block2 <- function(fit, model_label, keep_terms,pretty_map = c(Ve_z="$Ve$", sd_z="$sd$", logK_z="$log K$")
) {
  coefs <- tidy_fixed(fit, keep_terms)
  etas  <- eta_fixed(fit, keep_terms)
  r2s   <- r2(fit)
  out   <- merge(coefs, etas, by = "term", all = TRUE)
  out[, Term := pretty_map[term]]
  out[, `:=`(
    Model   = model_label,
    R2_m = unname(r2s$R2_marginal),
    R2_c = unname(r2s$R2_conditional)
  )]
  out
}
files <- list.files(run_root, recursive = TRUE, full.names = TRUE, pattern = "traitStats")
files <- files[grep("Hx1",files)]
# parse scenario metadata from path

summarize_century <- function(dt, t0, t1) {
  stopifnot(is.data.table(dt), "G" %in% names(dt))
  x <- dt[G >= t0 & G <= t1]
  if (nrow(x) == 0L) stop("No rows in the requested window.")

  # Per-generation derived metrics
  x[, `:=`(
    h2_LT   = Va_LT   / Vp_LT,
    IA_LT   = Va_LT   / (LT^2),
    h2_g    = Va_g    / Vp_g,
    IA_g    = Va_g    / (g^2),
    h2_hmax = Va_hmax / Vp_hmax,
    IA_hmax = Va_hmax / (hmax^2)
  )]

  # Helper: log-change between endpoints
  # v="LT"
  #plot(x$LT)
  get_r <- function(v) {
    v0 <- x[G == t0, get(v)][1]
    v1 <- x[G == t1, get(v)][1]
    if (is.na(v0) || is.na(v1) || v0 <= 0 || v1 <= 0) return(NA_real_)
    log(v1) - log(v0)
  }

  Tlen <- length(unique(x$G)) - 1L

  # Log and percent changes
  r_LT    <- get_r("LT")
  r_g     <- get_r("g")
  r_hmax  <- get_r("hmax")

  pct <- function(r) 100 * (exp(r) - 1)

  out <- data.table(
    trait     = c("LT","g","hmax"),
    r_total   = c(r_LT, r_g, r_hmax),
    r_per_gen = c(r_LT, r_g, r_hmax) / Tlen,
    pct_total = pct(c(r_LT, r_g, r_hmax)),
    pct_per_gen = pct(c(r_LT, r_g, r_hmax) / Tlen),
    h2_bar    = c(mean(x$h2_LT,   na.rm=TRUE),
                  mean(x$h2_g,    na.rm=TRUE),
                  mean(x$h2_hmax, na.rm=TRUE)),
    IA_bar    = c(mean(x$IA_LT,   na.rm=TRUE),
                  mean(x$IA_g,    na.rm=TRUE),
                  mean(x$IA_hmax, na.rm=TRUE)),
    VA_bar    = c(mean(x$Va_LT,   na.rm=TRUE),
                  mean(x$Va_g,    na.rm=TRUE),
                  mean(x$Va_hmax, na.rm=TRUE)),
    beta_bar  = c(mean(x$beta_LT,    na.rm=TRUE),
                  mean(x$beta_g,     na.rm=TRUE),
                  mean(x$beta_hmax,  na.rm=TRUE)),
    S_bar     = c(mean(x$S_LT,       na.rm=TRUE),
                  mean(x$S_g,        na.rm=TRUE),
                  mean(x$S_hmax,     na.rm=TRUE))
  )

  # Also return the per-generation table with derived quantities for plotting
  list(
    per_generation = x,
    window_summary = out,
    t0 = t0, t1 = t1, Tlen = Tlen
  )
}

# fp <- files[1]
# bname <- basename(fp)
run_summaries <- rbindlist(lapply(files,function(fp){
  #cat(fp,"\n")
  bname <<- basename(fp)

  pars <- data.table(
    Shmax = get_param("Shmax",bname),
    SLT   = get_param("SLT",bname),
    w      = get_param("w",bname),
    hlim  = get_param("hlim",bname),
    LTlim = get_param("LTlim",bname),
    f      = get_param("f",bname),
    Ve    = get_param("Ve",bname),
    sd    = get_param("sd",bname),
    Hx    = get_param("Hx",bname),
    HTarget      = get_param("H",bname),
    K      = get_param("K",bname)
  )

  dt <- fread(fp)
  if(max(dt$G)==11126){
    summary <- summarize_century(dt,t0=10000,t1=10125)
    cbind(pars,summary$window_summary)
  }else{
    return(NULL)
  }
}))
# Inspect
run_summaries[,theta_idx:=as.numeric(factor(paste(HTarget,Ve,sd,K,f,w,SLT,Shmax,Hx)))] ## each unique combination is a draw, each draw is simulated 5 times


run_summaries[,table(sd,Ve,K,Hx)]
# ---------- Trait labels as parsed expressions ----------
# Your dt has trait in {"LT","g","hmax"}.
dt <- copy(run_summaries)
dt[,CV_bar := sqrt(IA_bar)*100]
dt[,K := K/10000]
dt[, trait_lab := factor(trait,
                         levels = c("hmax","g","LT"),
                         labels = c("h[max]","g","L[T]"))]  # parsed by label_parsed

dt[sd==0.05 & Ve==0.05 & K==1 & HTarget==0.1270 & Hx==1,abs(mean(pct_total)),by=trait]

dt[sd==0.05 & Ve==0.05 & K==1& HTarget==0.1270,mean(CV_bar),by=trait]
dt[sd==0.05 & Ve==0.05 & K==1& HTarget==0.1270,mean(h2_bar),by=trait]
dt[sd==0.05 & Ve==0.05 & K==1 & HTarget==0.1270,mean(abs(pct_total)),by=trait][,mean(V1)]

dt[sd==0.05 & Ve==0.025 & K==1,mean(CV_bar),by=trait]
dt[sd==0.05 & Ve==0.025 & K==1,mean(h2_bar),by=trait]
dt[sd==0.05 & Ve==0.05 & K==1,abs(mean(pct_total)),by=trait][,mean(V1)]

dt[sd==0.05 & Ve==0.05 & K==1,mean(CV_bar)]
dt[sd==0.05 & Ve==0.05 & K==1,mean(h2_bar)]

#### ==== projected change ==== ####
run_root <- "./REEIM/scenarios/"
files <- list.files(run_root, recursive = TRUE, full.names = TRUE, pattern = "traitStats")
files <- files[grep("Hx1",files)]

run_summaries_proj <- rbindlist(lapply(files,function(fp){
  #cat(fp,"\n")
  bname <<- basename(fp)


  pars <- data.table(
    Shmax = get_param("Shmax",bname),
    SLT   = get_param("SLT",bname),
    w      = get_param("w",bname),
    hlim  = get_param("hlim",bname),
    LTlim = get_param("LTlim",bname),
    f      = get_param("f",bname),
    Ve    = get_param("Ve",bname),
    sd    = get_param("sd",bname),
    Hx    = get_param("Hx",bname),
    HTarget      = get_param("H",bname),
    K      = get_param("K",bname)
  )

  dt <- fread(fp)
  if(max(dt$G)==11126){
    summary <- summarize_century(dt,t0 = 10000,t1=10125)
    cbind(pars,summary$window_summary)
  }else{
    return(NULL)
  }
}))
# Inspect
run_summaries_proj[,theta_idx:=as.numeric(factor(paste(HTarget,Ve,sd,K,f,w,SLT,Shmax)))] ## each unique combination is a draw, each draw is simulated 2 times


#run_summaries_proj[,table(sd,Ve,K)]
# ---------- Trait labels as parsed expressions ----------
# Your dt has trait in {"LT","g","hmax"}.
dt <- copy(run_summaries_proj[HTarget==0.127])
dt[,CV_bar := sqrt(IA_bar)*100]
dt[,K := K/10000]
dt[, trait_lab := factor(trait,
                         levels = c("hmax","g","LT"),
                         labels = c("h[max]","g","L[T]"))]  # parsed by label_parsed
dt[sd==0.05 & Ve==0.05 & K==1 & HTarget==0.1270 & Hx==1,abs(mean(pct_total)),by=trait]
#dt[,unique(HTarget)]

dt[sd==0.05 & Ve==0.05 & K==1,mean(CV_bar),by=trait]
dt[sd==0.05 & Ve==0.05 & K==1,mean(h2_bar),by=trait]

dt[sd==0.05 & Ve==0.025 & K==1,mean(CV_bar),by=trait]
dt[sd==0.05 & Ve==0.025 & K==1,mean(h2_bar),by=trait]
dt[sd==0.05 & Ve==0.05 & K==1,abs(mean(pct_total)),by=trait][,mean(V1)]

dt[sd==0.05 & Ve==0.05 & K==1,mean(CV_bar)]
dt[sd==0.05 & Ve==0.05 & K==1,mean(h2_bar)]
#### ==== for SM ==== #####

traj_sum <- traj_base[, .(
  q05 = quantile(value, 0.05, na.rm = TRUE),
  q50 = quantile(value, 0.50, na.rm = TRUE),
  q95 = quantile(value, 0.95, na.rm = TRUE)
), by = .(Year, HTarget, Hx,Trait,sd,Ve,K)]
#traj_sum[,table(sd)]
traj_sum[,sd:=factor(sd)]
traj_sum[,Ve:=factor(Ve)]
traj_sum[,K:=factor(K)]

ylab_pct_total <- expression(Delta~"trait since onset, "~"%")
ylab_r_total   <- expression(Delta~"trait since onset, "~log~"scale")

# ---------- Generic panel builder (total change focus) ----------
make_panel <- function(dt, xvar, yvar, base_filter, xlab, ylab,
                       facet = trait_lab ~ .,
                       add_zero_line = FALSE) {
  dd <- copy(dt)[eval(base_filter)]
  stopifnot(nrow(dd) > 0)

  # Sample sizes for labels
  nlab <- dd[, .N, by = .(x = get(xvar), trait_lab, HTarget)]
  ymin <- dd[, .(ymin = min(get(yvar), na.rm = TRUE)), by = .(x = get(xvar), trait_lab, HTarget)]
  nlab <- nlab[ymin, on = .(x, trait_lab, HTarget)]
  nlab[, y := ymin - 0.02 * abs(ymin)]

  p <- ggplot(dd, aes(x = factor(get(xvar)), y = get(yvar))) +
    geom_violin(fill = "grey85", color = "grey40", trim = FALSE) +
    geom_boxplot(width = 0.18, outlier.shape = 16, outlier.size = 0.8, alpha = 0.6) +
    stat_summary(fun = median, geom = "point", size = 1.4, color = "black") +
    facet_grid(facet, scales = "free_y", labeller = labeller(.rows = label_parsed)) +
    labs(x = xlab, y = ylab) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      axis.text.x = element_text(angle = 0, vjust = 0.5)
    )

  if (add_zero_line) p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey40")

  # p + geom_text(data = nlab,
  #               aes(x = factor(x), y = y, label = paste0("n=", N)),
  #               size = 2.6, vjust = 1, inherit.aes = FALSE)
  p
}

# ---------- OFAT baselines ----------
base_sd <- quote(Ve == 0.05 & K == 1)   # vary sd
base_K  <- quote(Ve == 0.05 & sd == 0.05)   # vary K
base_Ve  <- quote(K == 1 & sd == 0.05)   # vary Ve

# ---------- Build panels: TOTAL change ----------
# Use pct_total (percent change over the century). If you'd rather log change, swap to yvar="r_total".

p_sd_pct <- make_panel(dt, xvar = "sd", yvar = "pct_total",
                       base_filter = base_sd,
                       xlab = "Mutational SD (sd)", ylab = ylab_pct_total,
                       add_zero_line = TRUE)

p_K_pct  <- make_panel(dt, xvar = "K",  yvar = "pct_total",
                       base_filter = base_K,
                       xlab = expression("Carrying capacity ("*K%*%10^5*")"), ylab = ylab_pct_total,
                       add_zero_line = TRUE)


p_Ve_pct  <- make_panel(dt, xvar = "Ve",  yvar = "pct_total",
                        base_filter = base_Ve,
                        xlab = expression("Environmental variance ("*V[e]*")"), ylab = ylab_pct_total,
                        add_zero_line = TRUE)

# Evolvability (mean over the century)
p_sd_IA <- make_panel(dt, xvar = "sd", yvar = "CV_bar",
                      base_filter = base_sd,
                      xlab = "Mutational SD (sd)", ylab = expression(bar(CV)[A]))

p_K_IA  <- make_panel(dt, xvar = "K",  yvar = "CV_bar",
                      base_filter = base_K,
                      xlab = expression("Carrying capacity ("*K%*%10^5*")"), ylab = expression(bar(CV)[A]))


p_Ve_IA  <- make_panel(dt, xvar = "Ve",  yvar = "CV_bar",
                       base_filter = base_Ve,
                       xlab = expression("Environmental variance ("*V[e]*")"), ylab = expression(bar(CV)[A]))



# Heritability (mean over the century)
p_sd_h2 <- make_panel(dt, xvar = "sd", yvar = "h2_bar",
                      base_filter = base_sd,
                      xlab = "Mutational SD (sd)", ylab = expression(bar(h)^2))

p_K_h2  <- make_panel(dt, xvar = "K",  yvar = "h2_bar",
                      base_filter = base_K,
                      xlab = expression("Carrying capacity ("*K%*%10^5*")"), ylab = expression(bar(h)^2))

p_Ve_h2 <- make_panel(dt, xvar = "Ve", yvar = "h2_bar",
                      base_filter = base_Ve,
                      xlab = expression("Environmental variance ("*V[e]*")"), ylab = expression(bar(h)^2))


# ---------- Titles and grid ----------

pdf("./Figures/Scenarios_summary.pdf",width = 12,height = 10)
(p_sd_pct+ggtitle("a") | p_K_pct+ggtitle("b") | p_Ve_pct+ggtitle("c")) /
  (p_sd_IA+ggtitle("d") | p_K_IA+ggtitle("e") | p_Ve_IA+ggtitle("f")) /
  (p_sd_h2+ggtitle("g") | p_K_h2+ggtitle("h") | p_Ve_h2+ggtitle("i"))
dev.off()
###################
### ETA SUMMARY ###
###################

# ---- 1) FIT MODELS (H_fac included as control variable) ----
dat_fit2 <- copy(run_summaries)
dat_fit2[,CV_bar := sqrt(IA_bar)*100]
dat_fit2[, `:=`(
  Ve_z   = as.numeric(scale(Ve)),
  sd_z   = as.numeric(scale(sd)),
  logK_z = as.numeric(scale(log(K))),
  H_fac  = factor(HTarget, levels = sort(unique(HTarget)))
)]


get_Rdem <- function(fit) {
  vc <- as.data.frame(VarCorr(fit))
  var_re <- vc[vc$grp == "theta_idx", "vcov"]
  var_res <- attr(VarCorr(fit), "sc")^2
  round(var_re / (var_re + var_res), 3)
}


# ---- 2) HELPER to extract Ve_z, sd_z, logK_z only (omit H terms) ----
tidy_fixed <- function(fit, keep_terms) {
  td <- as.data.table(tidy(fit, Parameters = "fixed", conf.int = TRUE))
  td <- td[term %in% keep_terms]
  setnames(td, c("estimate","conf.low","conf.high","std.error"),
           c("beta","ci_low","ci_high","se"))
  td[]
}

# ---- 3) HELPER to get partial eta^2 rows matching predictors ----
eta_fixed <- function(fit, keep_terms) {
  et <- as.data.table(eta_squared(fit, partial = TRUE, ci = 0.95,alternative = "two"))
  et <- et[Parameter %in% keep_terms]
  setnames(et, c("Parameter","Eta2_partial","CI_low","CI_high"),
           c("term","eta2p","eta2p_low","eta2p_high"))
  et[]
}

# ---- 4) Assemble table blocks for each model ----

assemble_block <- function(fit, model_label, keep_terms,pretty_map = c(Ve_z="$Ve$", sd_z="$sd$", logK_z="$log K$")
) {
  coefs <- tidy_fixed(fit, keep_terms)
  etas  <- eta_fixed(fit, keep_terms)
  r2s   <- r2(fit)
  out   <- merge(coefs, etas, by = "term", all = TRUE)
  out[, Term := pretty_map[term]]
  out[, `:=`(
    Model   = model_label,
    beta_fmt = sprintf("%.3f \\;[%.3f,\\; %.3f]", beta, ci_low, ci_high),
    eta_fmt  = sprintf("%.3f \\;[%.3f,\\; %.3f]", eta2p, eta2p_low, eta2p_high),
    R2_m = sprintf("%.3f", unname(r2s$R2_marginal)),
    R2_c = sprintf("%.3f", unname(r2s$R2_conditional))
  )]
  out[, .(Model, Term, `Std. $\\beta$ [95\\% CI]` = beta_fmt,
          `Partial $\\eta^2$ [95\\% CI]` = eta_fmt,
          `$R^2_m$` = R2_m, `$R^2_c$` = R2_c)]
}

assemble_block2 <- function(fit, model_label, keep_terms,pretty_map = c(Ve_z="$Ve$", sd_z="$sd$", logK_z="$log K$")
) {
  coefs <- tidy_fixed(fit, keep_terms)
  etas  <- eta_fixed(fit, keep_terms)
  r2s   <- r2(fit)
  out   <- merge(coefs, etas, by = "term", all = TRUE)
  out[, Term := pretty_map[term]]
  out[, `:=`(
    Model   = model_label,
    R2_m = unname(r2s$R2_marginal),
    R2_c = unname(r2s$R2_conditional)
  )]
  out
}


tab_traits <- rbindlist(lapply(c("g","LT","hmax"),function(tr){
  print(tr)
  fit_h2 <- lmer(h2_bar    ~ Ve_z + sd_z + logK_z + H_fac + (1|theta_idx), data = dat_fit2[trait==tr])
  fit_CV <- lmer(CV_bar    ~ Ve_z + sd_z + logK_z + H_fac + (1|theta_idx), data = dat_fit2[trait==tr])
  fit_r  <- lmer(r_total ~ h2_bar + CV_bar  + H_fac + (1|theta_idx), data = dat_fit2[trait==tr])


  keep <- c("Ve_z","sd_z","logK_z")
  tab_h2 <- assemble_block(fit_h2, model_label = "Mean heritability ($\\overline{h^2}$)", keep,pretty_map = c(Ve_z="$Ve$", sd_z="$sd$", logK_z="$log K$"))
  tab_CV <- assemble_block(fit_CV, model_label = "Mean evolvability ($\\overline{I_A}$)", keep,pretty_map = c(Ve_z="$Ve$", sd_z="$sd$", logK_z="$log K$"))
  keep_terms <- keep <- c("h2_bar", "CV_bar")
  tab_r <- assemble_block(fit_r,"Log response. ($r_T$)", keep,pretty_map = c(h2_bar="$\\overline{h^2}$", CV_bar="$\\overline{CV_A}$"))
  tab_all <- rbindlist(list(tab_h2, tab_CV, tab_r), use.names = TRUE, fill = TRUE)

  # print R2 only once per model block
  tab_all[, `$R^2_m$` := ifelse(duplicated(Model), "", `$R^2_m$`)]
  tab_all[, `$R^2_c$` := ifelse(duplicated(Model), "", `$R^2_c$`)]
  tab_all[,trait:=tr]
  return(tab_all)

}))

tab_traits <- tab_traits[,.( Model,Term,trait,`Std. $\\beta$ [95\\% CI]`,`Partial $\\eta^2$ [95\\% CI]`,`$R^2_m$`,`$R^2_c$`)]
tab_traits[trait=="g",trait:="$g$"]
tab_traits[trait=="LT",trait:="$L_T$"]
tab_traits[trait=="hmax",trait:="$h_{max}$"]

#c(g="g",hmax="h[max]",LT="L[T]")
# ---- 5) LaTeX table (H included but not interpreted) ----
cap <- paste(
  "Standardized Parameters of environmental variance (Ve), mutational variance (sd),",
  "and carrying capacity (log K) on mean heritability ($\\overline{h^2}$), mean coefficient of additive genetic variation",
  "($\\overline{CV_A}$), and per-generation log response ($R/T$).",
  "Harvest level ($H_\\mathrm{target}$) was included as a fixed factor in all models to control",
  "for mean differences among fishing intensities, but its coefficients were not interpreted.",
  "Entries are standardized coefficients (predictors z-scored) with 95\\% CIs, partial $\\eta^2$ with 95\\% CIs,",
  "and marginal/conditional $R^2$."
)

kbl <- kable(tab_traits,
             format = "latex", booktabs = TRUE, escape = FALSE,
             col.names = c("Model","Predictor","Trait","Std. $\\beta$ [95\\% CI]",
                           "Partial $\\eta^2$ [95\\% CI]","$R^2_m$","$R^2_c$"),
             caption = cap, label = "Parametersizes_main") |>
  kable_styling(latex_options = c("hold_position")) |>
  collapse_rows(columns = 1, latex_hline = "none") |>
  add_footnote(
    label = "$H_\\mathrm{target}$ (low, medium, high) included as a categorical control; coefficients and interactions not shown.",
    notation = "symbol",
    escape = FALSE
  )
#kbl
save_kable(kbl, file = "./Tables/Table_Parametersizes_main.tex")


assemble_block_numeric <- function(fit, model_label, keep_terms,pretty_map = c(Ve_z="$Ve$", sd_z="$sd$", logK_z="$log K$")
) {
  coefs <- tidy_fixed(fit, keep_terms)
  etas  <- eta_fixed(fit, keep_terms)
  r2s   <- r2(fit)
  out   <- merge(coefs, etas, by = "term", all = TRUE)
  out[, Term := pretty_map[term]]
  out[, `:=`(
    Model   = model_label,
    R2_m = unname(r2s$R2_marginal),
    R2_c = unname(r2s$R2_conditional)
  )]
  out
}
#####


tab_traits_numeric <- rbindlist(lapply(c("g","LT","hmax"),function(tr){
  print(tr)
  fit_h2 <- lmer(h2_bar    ~ Ve_z + sd_z + logK_z + H_fac + (1|theta_idx), data = dat_fit2[trait==tr])
  fit_CV <- lmer(CV_bar    ~ Ve_z + sd_z + logK_z + H_fac + (1|theta_idx), data = dat_fit2[trait==tr])
  fit_r  <- lmer(r_total ~ h2_bar + CV_bar + logK_z + H_fac + (1|theta_idx), data = dat_fit2[trait==tr])


  keep <- c("Ve_z","sd_z","logK_z")
  tab_h2 <- assemble_block_numeric(fit_h2, model_label = '"a)"~ h^2 %~% V[e] + sd + log(K)', keep,pretty_map = c(Ve_z="V[e]", sd_z="sd", logK_z="log(K)"))
  tab_CV <- assemble_block_numeric(fit_CV, model_label = '"b)"~ CV[A] %~% V[e] + sd + log(K)', keep,pretty_map = c(Ve_z="V[e]", sd_z="sd", logK_z="log(K)"))
  keep_terms <- keep <- c("h2_bar", "CV_bar","logK_z")
  tab_r <- assemble_block2(fit_r,model_label ='"c)"~ r[T] %~% h^2 + CV[A] + log(K)', keep,pretty_map = c(h2_bar="h^2", CV_bar="CV[A]", logK_z="log(K)"))
  tab_all <- rbindlist(list(tab_h2, tab_CV, tab_r), use.names = TRUE, fill = TRUE)

  # print R2 only once per model block
  tab_all[,trait:=tr]
  return(tab_all)

}))

dt <- tab_traits_numeric[,.(Model,Term,trait,eta2p,eta2p_low,  eta2p_high)]

dt[, Term := factor(Term, levels = c("V[e]","sd","log(K)","h^2","CV[A]"))]

pd <- position_dodge(width = 0.7)

p_eta_summary <- ggplot(dt, aes(x = trait, y = eta2p, fill = Term)) +
  geom_col(position = pd, width = 0.6, color = "grey25") +
  geom_errorbar(aes(ymin = eta2p_low, ymax = eta2p_high),
                position = pd, width = 0.15, color = "grey20") +
  facet_grid(. ~ Model, labeller = labeller(.cols = label_parsed)) +
  scale_fill_manual(values = alpha(pal_jco()(10), 1),
                    labels = function(x) parse(text = x),
                    name = "Predictor") +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0,1),
                     name = expression(eta[partial]^2)) +
  scale_x_discrete(labels = function(x) parse(text = c(g="g", hmax="h[max]", LT="L[T]")[x]),
                   name = NULL) +
  theme_bw() +
  theme(aspect.ratio = 1,
        strip.text.x = element_text(hjust = 0),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom")

pdf("./Figures/eta_sum.pdf",width = 8,height = 3.5)
p_eta_summary
dev.off()


dt_r <- copy(dt)

# Ensure stable legend order
dt_r[, Term := factor(Term, levels = c("V[e]","sd","log(K)","h^2","CV[A]"))]

# Panel-specific ranking: largest eta2p first per (Model, trait)
dt_r[, rank_within := -rank(eta2p, ties.method = "first"),
     by = .(Model, trait)]
# Build a panel-specific x var combining trait and ranked Term
dt_r[, x_panel := interaction(trait, Term, drop = TRUE)]
# To display grouped bars by Term in trait order, we keep x as trait, but order the legend by Term; for
# within-panel ranking on the x-axis, switch to this:
dt_r[, Term_panel := factor(Term, levels = unique(Term[order(Model, trait, rank_within)]))]


pd <- position_dodge(width = 0.7)

p_eta_ranked <- ggplot(dt_r,
                       aes(x = trait, y = eta2p, fill = Term_panel)) +
  geom_col(position = pd, width = 0.6, color = "grey25") +
  geom_errorbar(aes(ymin = eta2p_low, ymax = eta2p_high),
                position = pd, width = 0.15, color = "grey20") +
  facet_grid(. ~ Model, labeller = labeller(.cols = label_parsed)) +
  scale_fill_manual(
    values = alpha(pal_jco()(10), 1),
    labels = function(x) parse(text = as.character(x)),
    name = "Predictor"
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0,1),
                     name = expression(eta[partial]^2)) +
  scale_x_discrete(labels = function(x) parse(text = c(g="g", hmax="h[max]", LT="L[T]")[x]),
                   name = NULL) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    strip.text.x = element_text(hjust = 0),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom"
  )


pdf("./Figures/eta_sum.pdf",width = 8,height = 3.5)
p_eta_ranked
dev.off()
