# ---- PACKAGES ----
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

# Modeling
library(lme4)        # mixed-effects
library(broom.mixed) # tidy summaries (optional)


run_root <- "./REEIM/scenarios/"      # folder that contains scenario subfolders
t0 <- 10000; t1 <- 10125                      # century window
trait_names <- c("LT","g","hmax")             # traits to summarize
cores <- min(10, detectCores())

# list all trait_stats files

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
                  mean(x$Va_hmax, na.rm=TRUE))
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
    SHmax = get_param("Shmax",bname),
    SLT   = get_param("SLT",bname),
    w      = get_param("w",bname),
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
run_summaries[,theta_idx:=as.numeric(factor(paste(HTarget,Ve,sd,K,f,w,SLT,SHmax)))] ## each unique combination is a draw, each draw is simulated 5 times
#run_summaries[,table(sd)]


# ---------- Trait labels as parsed expressions ----------
# Your dt has trait in {"LT","g","hmax"}.
dt <- copy(run_summaries)
dt[,CV_bar := sqrt(IA_bar)*100]
dt[,K := K/10000]
dt[, trait_lab := factor(trait,
                         levels = c("hmax","g","LT"),
                         labels = c("h[max]","g","L[T]"))]  # parsed by label_parsed



# ---------- Y labels ----------
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
dt[sd==0.05 & Ve==0.05 & K==1,mean(CV_bar)]
dt[sd==0.05 & Ve==0.05 & K==1,mean(h2_bar)]

# ---------- Titles and grid ----------


row1 <-  grid.arrange(p_sd_pct+ggtitle("a"), p_K_pct+ggtitle("b"), p_Ve_pct+ggtitle("c"),ncol = 3)

row2 <- grid.arrange(p_sd_IA+ggtitle("d"), p_K_IA+ggtitle("e"), p_Ve_IA+ggtitle("f"),ncol = 3)

row3 <- grid.arrange(p_sd_h2+ggtitle("g"), p_K_h2+ggtitle("h"), p_Ve_h2+ggtitle("i"),ncol = 3)

pdf("./Figures/Scenarios_summary.pdf",width = 12,height = 10)
grid.arrange(row1, row2, row3, ncol = 1)
dev.off()


# ---- 1) FIT MODELS (H_fac included as control variable) ----
dat_fit2 <- copy(run_summaries)
dat_fit2[,CV_bar := sqrt(IA_bar)*100]
dat_fit2[, `:=`(
  Ve_z   = as.numeric(scale(Ve)),
  sd_z   = as.numeric(scale(sd)),
  logK_z = as.numeric(scale(log(K))),
  H_fac  = factor(HTarget, levels = sort(unique(HTarget)))
)]

# tmp <- dat_fit2[,table(theta_idx)]
# keep_idx <- as.numeric(names(tmp)[tmp>3])
# dat_fit2 <- dat_fit2[theta_idx %in%keep_idx]

get_Rdem <- function(fit) {
  vc <- as.data.frame(VarCorr(fit))
  var_re <- vc[vc$grp == "theta_idx", "vcov"]
  var_res <- attr(VarCorr(fit), "sc")^2
  round(var_re / (var_re + var_res), 3)
}

### these are all zero!
# R_dem_h2 <- get_Rdem(fit_h2)
# R_dem_IA <- get_Rdem(fit_IA)
# R_dem_R  <- get_Rdem(fit_r)

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
#eta
# ---- 4) Assemble table blocks for each model ----
#fit <- fit_r
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
#tr = "g"
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
  #summary(fit_r)
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
save_kable(kbl, file = "Table_Parametersizes_main.tex")


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

#####
tab_traits_numeric <- rbindlist(lapply(c("g","LT","hmax"),function(tr){
  print(tr)
  fit_h2 <- lmer(h2_bar    ~ Ve_z + sd_z + logK_z + H_fac + (1|theta_idx), data = dat_fit2[trait==tr])
  fit_CV <- lmer(CV_bar    ~ Ve_z + sd_z + logK_z + H_fac + (1|theta_idx), data = dat_fit2[trait==tr])
  fit_r  <- lmer(r_total ~ h2_bar + CV_bar  + H_fac + (1|theta_idx), data = dat_fit2[trait==tr])
  
  
  keep <- c("Ve_z","sd_z","logK_z")
  tab_h2 <- assemble_block_numeric(fit_h2, model_label = '"a)"~ h^2 %~% V[e] + sd + log(K)', keep,pretty_map = c(Ve_z="V[e]", sd_z="sd", logK_z="log(K)"))
  tab_CV <- assemble_block_numeric(fit_CV, model_label = '"b)"~ CV[A] %~% V[e] + sd + log(K)', keep,pretty_map = c(Ve_z="V[e]", sd_z="sd", logK_z="log(K)"))
  keep_terms <- keep <- c("h2_bar", "CV_bar")
  tab_r <- assemble_block2(fit_r,model_label ='"c)"~ r[T] %~% h^2 + CV[A] + log(K)', keep,pretty_map = c(h2_bar="h^2", CV_bar="CV[A]"))
  #summary(fit_r)
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

#### evolutionary change ####
st all trait_stats files

files <- list.files(run_root, recursive = TRUE, full.names = TRUE, pattern = "H*trait_stats")
files <- files[grep("Hx1",files)]
files[1]
fi <- 1

gens <- c(seq(9900,10325,by=5))
get_param <- function(key,bname) {
  as.numeric(strsplit(strsplit(bname,key)[[1]][2],"_",fixed = TRUE)[[1]][1])
}

out <- rbindlist(mclapply(1:length(files),function(fi){
  print(fi)
  fp <- files[fi]
  
  bname <- basename(fp)
  
  pars <- data.table(
    Ve    = get_param("Ve_",bname),
    sd    = get_param("sd_",bname),
    HTarget      = get_param("H_",bname),
    K      = get_param("K_",bname)
  )
  dt <- fread(fp)
  out <- cbind(pars,dt[G%in%gens,.(G,LT,g,hmax,h,median_size)],theta_idx=fi)
  
},mc.cores = 1))   



grid <- as.matrix(unique(out[,.(K,sd,Ve,HTarget,theta_idx)]))
mean_LT <- 169.2
mean_hmax <- 62
mean_h <- 38
mean_g <- 0.24
#€mean_size <- out[Year==2025,mean(size)]


out2 <- rbindlist(apply(grid,1,function(x){
  tmp <- out[K==x[1] & sd==x[2] & Ve==x[3] & HTarget==x[4] & theta_idx==x[5]]
  
  delta <- mean(tmp[G==10125,LT])-mean_LT
  tmp[,LT:=LT-delta]
  
  delta <- mean(tmp[G==10125,g])-mean_g
  tmp[,g:=g-delta]
  
  delta <- mean(tmp[G==10125,hmax])-mean_hmax
  tmp[,hmax:=hmax-delta]
  
  delta <- mean(tmp[G==10125,h])-mean_h
  tmp[,h:=h-delta]
  
}))

out2[,Year:=G-10000+1900]

dt2 <- melt(out2,measure.vars = c("h","g","LT"))

dt2[variable=="h",Trait:="Annual juvenile growth increment (mm)"]
dt2[variable=="LT",Trait:="Size at maturation (mm)"]
dt2[variable=="g",Trait:="Parental investment"]

trait_evolution_data <- copy(dt2)

dt <- trait_evolution_data[K==10000 & Ve==0.05 & sd==0.05 & Year>=1800 & Year<=2225] ## default settings
dt <- dt[!is.na(G),.(value=mean(value)),by=.(G,HTarget,Trait,Year)]

ggplot(dt,aes(Year,value,col=factor(HTarget)))+
  geom_line() +
  facet_wrap(.~Trait,scales="free_y",nrow=3)+
  theme_bw(base_size = 15) +
  scale_color_manual(
    values = alpha(pal_jco()(10)[c(1,2,4)], 1),
    labels = function(x) parse(text = x),
    name="Harvest rate (% biomass per year)"
  ) +
  theme(
    aspect.ratio = 0.75,
    strip.text.x = element_text(hjust = 0),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom"
  ) +
  geom_vline(xintercept = 1900,linetype=2)+
  geom_vline(xintercept = 2025,linetype=2) +
  xlab("Year")+
  ylab("Value")



dt3 <- melt(out2,measure.vars = c("h","g","LT"))

dt3[variable=="h",Trait:="Annual juvenile growth increment (mm)"]
dt3[variable=="LT",Trait:="Size at maturation (mm)"]
dt3[variable=="g",Trait:="Parental investment"]

trait_evol_full <- dt3
saveRDS(dt3,"trait_evol_full.rds")
trait_evol_full <- readRDS("trait_evol_full.rds")
dt <- trait_evol_full[K==10000 & Ve==0.05 & sd==0.05 & Year>=1800 & Year<=2225 & HTarget==0.032] ## default settings
dt2 <- dt[,.(value=mean(value)),by=.(Year,Trait)]

ggplot(dt2,aes(Year,value))+
  geom_line(linewidth = 0.15,alpha=0.5,col="grey30") +
  geom_line(data=dt2,linewidth = 2,alpha=1,aes(Year,value),inherit.aes = FALSE) +
  facet_wrap(.~Trait,scales="free_y",nrow=3)+
  
  theme_bw(base_size = 20) +
  scale_color_manual(
    values = alpha(pal_jco()(10)[c(1,2,4)], 1),
    labels = function(x) parse(text = x),
    name="Harvest rate (% biomass per year)"
  ) +
  theme(
    aspect.ratio = 0.75,
    strip.text.x = element_text(hjust = 0),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom"
  ) +
  geom_vline(xintercept = 1900,linetype=2)+
  geom_vline(xintercept = 2025,linetype=2) +
  xlab("Year")+
  ylab("Value")


