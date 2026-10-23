library(data.table); library(ggplot2); t0 <- 10000; t1 <- 10125
base_folder <- "/private/tmp/claude-539526166/-Users-petrikem-gitlab-REEIM/f7799891-e671-4825-bad9-a218acb19f48/scratchpad/final/res_plaice/"
APOT        <- 10125   
REEIM_res <- readRDS(file.path(base_folder, "REEIM_res.rds"))
acc   <- REEIM_res$particles[[5]][which(keep)]

run_root <- "/private/tmp/claude-539526166/-Users-petrikem-gitlab-REEIM/f7799891-e671-4825-bad9-a218acb19f48/scratchpad/final/res_plaice/GEN_5/"
files <- list.files(run_root, recursive = TRUE, full.names = TRUE,pattern = "traitStats")


## ------------------------------------------------------------
## 1) Century summarizer (signed r_T; heritability/evolvability)
## ------------------------------------------------------------
summarize_century <- function(dt, t0, t1) {
  stopifnot(is.data.table(dt), "G" %in% names(dt))
  x <- dt[G >= t0 & G <= t1]
  if (nrow(x) == 0L) stop("No rows in the requested window.")
  
  # Derived per-generation quantities
  x[, `:=`(
    h2_LT   = Va_LT   / Vp_LT,
    IA_LT   = Va_LT   / (LT^2),
    h2_g    = Va_g    / Vp_g,
    IA_g    = Va_g    / (g^2),
    h2_hmax = Va_hmax / Vp_hmax,
    IA_hmax = Va_hmax / (hmax^2)
  )]
  
  # Helper: log change between endpoints
  get_r <- function(v) {
    v0 <- x[G == t0, get(v)][1]
    v1 <- x[G == t1, get(v)][1]
    if (is.na(v0) || is.na(v1) || v0 <= 0 || v1 <= 0) return(NA_real_)
    log(v1) - log(v0)
  }
  
  Tlen <- t1 - t0          # explicit century length in generations
  r_LT    <- get_r("LT")
  r_g     <- get_r("g")
  r_hmax  <- get_r("hmax")
  
  pct <- function(r) 100 * (exp(r) - 1)
  
  out <- data.table(
    trait       = c("LT","g","hmax"),
    r_total     = c(r_LT, r_g, r_hmax),
    r_per_gen   = c(r_LT, r_g, r_hmax) / Tlen,
    pct_total   = pct(c(r_LT, r_g, r_hmax)),
    pct_per_gen = pct(c(r_LT, r_g, r_hmax) / Tlen),
    h2_bar      = c(mean(x$h2_LT,   na.rm=TRUE),
                    mean(x$h2_g,    na.rm=TRUE),
                    mean(x$h2_hmax, na.rm=TRUE)),
    IA_bar      = c(mean(x$IA_LT,   na.rm=TRUE),
                    mean(x$IA_g,    na.rm=TRUE),
                    mean(x$IA_hmax, na.rm=TRUE)),
    VA_bar      = c(mean(x$Va_LT,   na.rm=TRUE),
                    mean(x$Va_g,    na.rm=TRUE),
                    mean(x$Va_hmax, na.rm=TRUE))
  )
  list(per_generation = x, window_summary = out, Tlen = Tlen)
}

## ------------------------------------------------------------
## 2) Run summaries for all accepted trait files (5th gen)
## ------------------------------------------------------------

run_summaries_base <- rbindlist(lapply(seq_len(nrow(acc)), function(i){
 
  
  dt <- fread(files[grep(acc$file[i],files)], showProgress = FALSE)
  if (!("G" %in% names(dt)) || max(dt$G, na.rm=TRUE) < t1) return(NULL)
  s <- summarize_century(dt, t0, t1)$window_summary
  # attach particle-level params (H, weights if present)
  s[, H := acc$H[i]]
  if ("weights" %in% names(acc)) s[, weight := acc$weights[i]]
  s[, file := acc$file[i]]
  s
}), fill = TRUE)

# Tidy trait names for plotting: hmax -> h, LT -> L[T] as parsed math
run_summaries_base[, trait_plot := fifelse(trait == "hmax","A~~growth~(h[max])",fifelse(trait == "LT","C~~length~at~maturation~(L[T])","B~~reproductive~investment~(g)"))]
run_summaries_base[, trait_plot := factor(trait_plot, levels = c("A~~growth~(h[max])","B~~reproductive~investment~(g)","C~~length~at~maturation~(L[T])"))]


## ------------------------------------------------------------
## 3) Correlations and simple slopes (per trait)
## ------------------------------------------------------------
# Use signed r_per_gen; if you prefer percent per gen, swap to pct_per_gen/100
dat_corr <- run_summaries_base[!is.na(H) & !is.na(pct_total)]

p <- ggplot(dat_corr,
            aes(x = H, y = pct_total, group = trait_plot)) +
  geom_point(alpha = 0.6, size = 1.8,col="grey50") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8,col="black") +
  facet_grid(. ~ trait_plot, labeller = label_parsed) +
  scale_x_continuous("harvest intensity H (fraction of biomass removed per year)", limits = c(0, NA)) +
  scale_y_continuous(expression(r[T]~", %"),
                     labels = function(v) sprintf("%.1f", v)) +
  theme_bw(base_size = 11) +
  theme(
    aspect.ratio = 1,
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )


# corr_tbl <- dat_corr[, {
#   # Pearson & Spearman
#   p_cor <- suppressWarnings(cor.test(H, pct_total, method = "pearson"))
#   s_cor <- suppressWarnings(cor.test(H, pct_total, method = "spearman"))
#   # Simple LM slope with 95% CI
#   fit <- lm(r_per_gen ~ H)
#   cf  <- coef(summary(fit))["H", ]
#   ci  <- confint(fit)["H", ]
#   data.table(
#     n               = .N,
#     pearson_r       = unname(p_cor$estimate),
#     pearson_p       = p_cor$p.value,
#     spearman_rho    = unname(s_cor$estimate),
#     spearman_p      = s_cor$p.value,
#     slope_H         = unname(cf["Estimate"]),
#     slope_se        = unname(cf["Std. Error"]),
#     slope_ci_lo     = ci[1],
#     slope_ci_hi     = ci[2]
#   )
# }, by = trait_plot]
# 
# # Neat LaTeX table
# corr_tbl_disp <- copy(corr_tbl)
# corr_tbl_disp[, `Slope (H → r_T/T) [95% CI]` :=
#                 sprintf("%.4f [%.4f, %.4f]", slope_H, slope_ci_lo, slope_ci_hi)]
# corr_tbl_disp[, `Pearson r (p)`   := sprintf("%.3f (%.3g)", pearson_r,  pearson_p)]
# corr_tbl_disp[, `Spearman ρ (p)`  := sprintf("%.3f (%.3g)", spearman_rho, spearman_p)]
# corr_tbl_disp <- corr_tbl_disp[, .(Trait = trait_plot, n,
#                                    `Pearson r (p)`, `Spearman ρ (p)`,
#                                    `Slope (H → r_T/T) [95% CI]`)]
# 
# kbl <- kable(corr_tbl_disp, format = "latex", booktabs = TRUE, escape = FALSE,
#              caption = "Correlation between harvest rate ($H$) and per-generation log change ($r_T/T$) among accepted ABC particles (5th generation), by trait.",
#              col.names = c("Trait","n","Pearson $r$ (p)","Spearman $\\rho$ (p)","Slope $\\beta_H$ [95\\% CI]"),
#              label = "tab:cor_H_rT") |>
#   kable_styling(latex_options = "hold_position")
# save_kable(kbl, "Table_Correlation_H_rT.tex")

## ------------------------------------------------------------
## 4) Plot: H vs r_T/T (signed), by trait
## ------------------------------------------------------------
w_tbl <- dat_corr[
  , {
    fitw <- lm(r_total ~ H, weights = weight)
    ci   <- confint(fitw)["H", ]
    est  <- coef(summary(fitw))["H", "Estimate"]
    # Weighted Pearson r via cov_w / sd_w decomposition
    w    <- if ("weight" %in% names(.SD)) weight else rep(1, .N)
    w    <- w / sum(w)
    Hp   <- H - sum(w*H); Rp <- r_per_gen - sum(w*r_per_gen)
    r_w  <- sum(w*Hp*Rp) / sqrt(sum(w*Hp^2) * sum(w*Rp^2))
    data.table(slope_H = est, slope_ci_lo = ci[1], slope_ci_hi = ci[2], pearson_r = r_w)
  },
  by = trait_plot
]

corr_annot_w <- w_tbl[
  , .(trait_plot,
      label = sprintf("italic(r)==%.2f~~beta[H]==%.3f",
                      pearson_r, slope_H))
]

p_annot_w <- p +
  geom_text(
    data = corr_annot_w,
    aes(x = 0.003, y = Inf, label = label),
    parse = TRUE, hjust = 0, vjust = 1.5, size = 3.5
  )

ggsave("./compare_old_new/figures_dardel_plaice/Evol_resp_vs_H_plaice.pdf", p_annot_w, height = 4, width = 12); ggsave("./compare_old_new/figures_dardel_plaice/Evol_resp_vs_H_plaice.png", p_annot_w, height = 4, width = 12, dpi = 200)


