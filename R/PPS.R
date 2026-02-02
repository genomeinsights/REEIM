library(ggplot2)
library(ggsci)
library(data.table)
run_root <- "./REEIM/scenarios3/"
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

#### === violin plots of scenarios === ###


run_root <- "./REEIM/scenarios3/"      # folder that contains scenario subfolders
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
#files <- files[grep("Hx1",files)]
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
  get_r <- function(v) {
    v0 <- x[G == t0, get(v)][1]
    v1 <- x[G == t1, get(v)][1]
    if (is.na(v0) || is.na(v1) || v0 <= 0 || v1 <= 0) return(NA_real_)
    log(v1) - log(v0)
  }
  
  Tlen <- length(unique(x$G)) - 1L  # should be 100 if generations are contiguous
  
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
  cat(fp,"\n")
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
    HTarget      = get_param("H",bname),
    K      = get_param("K",bname)
  )
  
  dt <- fread(fp)
  if(max(dt$G)==11126){
    summary <- summarize_century(dt,t0,t1)
    cbind(pars,summary$window_summary)  
  }else{
    return(NULL)
  }
}))
# Inspect
run_summaries[,theta_idx:=as.numeric(factor(paste(HTarget,Ve,sd,K,f,w,SLT,Shmax)))] ## each unique combination is a draw, each draw is simulated 5 times


run_summaries[,table(sd,Ve,K)]
# ---------- Trait labels as parsed expressions ----------
# Your dt has trait in {"LT","g","hmax"}.
dt <- copy(run_summaries)
dt[,CV_bar := sqrt(IA_bar)*100]
dt[,K := K/10000]
dt[, trait_lab := factor(trait,
                         levels = c("hmax","g","LT"),
                         labels = c("h[max]","g","L[T]"))]  # parsed by label_parsed




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

ylab_pct_total <- expression(Delta~"trait since onset, "~"%")
ylab_r_total   <- expression(Delta~"trait since onset, "~log~"scale")

# ---------- Generic panel builder (total change focus) ----------
make_panel <- function(dt, xvar, yvar, base_filter, xlab, ylab,
                       facet = trait_lab ~ paste("H =",HTarget),
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
      aspect.ratio = 0.75,
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

dt[sd==0.05 & Ve==0.05 & K==1,mean(CV_bar),by=trait]
dt[sd==0.05 & Ve==0.05 & K==1,mean(h2_bar),by=trait]
dt[sd==0.05 & Ve==0.05 & K==1,abs(mean(pct_total)),by=trait][,mean(V1)]

dt[sd==0.05 & Ve==0.025 & K==1,mean(CV_bar),by=trait]
dt[sd==0.05 & Ve==0.025 & K==1,mean(h2_bar),by=trait]
dt[sd==0.05 & Ve==0.05 & K==1,abs(mean(pct_total)),by=trait][,mean(V1)]

dt[sd==0.05 & Ve==0.05 & K==1,mean(CV_bar)]
dt[sd==0.05 & Ve==0.05 & K==1,mean(h2_bar)]

# ---------- Titles and grid ----------


row1 <-  grid.arrange(p_sd_pct+ggtitle("a"), p_K_pct+ggtitle("b"), p_Ve_pct+ggtitle("c"),ncol = 3)

row2 <- grid.arrange(p_sd_IA+ggtitle("d"), p_K_IA+ggtitle("e"), p_Ve_IA+ggtitle("f"),ncol = 3)

row3 <- grid.arrange(p_sd_h2+ggtitle("g"), p_K_h2+ggtitle("h"), p_Ve_h2+ggtitle("i"),ncol = 3)

pdf("./Figures/Scenarios_summary.pdf",width = 12,height = 10)
grid.arrange(row1, row2, row3, ncol = 1)
dev.off()

