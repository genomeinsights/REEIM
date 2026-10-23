## Perch PPS from the final posterior (GEN_6), Dardel jobs 26792544-6. Adapted from R/PPS.R.
## Args: folder with the scenario subfolders (traitStats only), output folder.
suppressMessages({library(data.table); library(ggplot2); library(lme4); library(effectsize); library(performance); library(patchwork); library(ggsci); library(scales)})
a <- commandArgs(TRUE); run_root <- a[1]; out <- a[2]
gp <- function(key, b) as.numeric(strsplit(strsplit(b, key)[[1]][2], "_", fixed = TRUE)[[1]][1])
files <- list.files(run_root, recursive = TRUE, full.names = TRUE, pattern = "^traitStats_")
APOT <- 10125; G_end <- 11126
tlab <- c(h = "A~~annual~growth~increment~(h)", LT = "B~~length~at~maturation~(L[T])", g = "C~~reproductive~investment~(g)")
## ---- read all runs once ----
runs <- rbindlist(lapply(files, function(fp) { b <- basename(fp); dt <- fread(fp)
  data.table(file = b, Ve = gp("Ve", b), sd = gp("sd", b), K = gp("K", b), HTarget = gp("_H", b), Hx = gp("Hx", b),
             w = gp("_w", b), f = gp("_f", b), SLT = gp("SLT", b), Shmax = gp("Shmax", b), dt = list(dt)) }))
runs[, complete := vapply(dt, function(d) max(d$G) == G_end, logical(1))]
cat("runs:", nrow(runs), " complete:", sum(runs$complete), "\n"); runs <- runs[complete == TRUE]
runs[, theta := as.integer(factor(paste(HTarget, Ve, sd, K, Hx, w, f, SLT, Shmax)))]
base <- quote(Ve == 0.05 & sd == 0.05 & K == 10000)
## ---- 1) trajectories (Fig. 3 / SI), recentred so the 2025 mean equals the empirical values ----
ref <- c(LT = 169.2, hmax = 62, h = 38, g = 0.24); gens <- seq(9900, 10325, by = 5)
traj <- runs[eval(base), { d <- dt[[1]][G %in% gens, .(G, LT, g, h)]; d }, by = .(HTarget, Hx, theta, file)]
traj[, `:=`(LT = LT - (mean(LT[G == APOT]) - ref["LT"]), g = g - (mean(g[G == APOT]) - ref["g"]), h = h - (mean(h[G == APOT]) - ref["h"])), by = .(HTarget, Hx)]
traj[, Year := G - 10000 + 1900]
tl <- melt(traj, id.vars = c("HTarget", "Hx", "theta", "file", "Year"), measure.vars = c("h", "LT", "g"))
ts <- tl[, .(q05 = quantile(value, .05), q50 = median(value), q95 = quantile(value, .95)), by = .(Year, HTarget, Hx, variable)]
ts[, Trait := factor(tlab[as.character(variable)], levels = tlab)]
ts[, HT := factor(sprintf("H == %.2f", HTarget), levels = sprintf("H == %.2f", sort(unique(HTarget))))]
fwrite(ts, file.path(out, "pps_trajectory_quantiles.csv"))
trajplot <- function(d, facet) ggplot(d, aes(Year)) +
  geom_ribbon(data = d[Hx == 1], aes(ymin = q05, ymax = q95), alpha = .25, fill = "firebrick") +
  geom_ribbon(data = d[Hx == 0 & Year >= 2025], aes(ymin = q05, ymax = q95), alpha = .25, fill = "steelblue") +
  geom_line(data = d[Hx == 1], aes(y = q50), linewidth = 1, colour = "firebrick") +
  geom_line(data = d[Hx == 0 & Year >= 2025], aes(y = q50), linewidth = .75, colour = "steelblue") +
  geom_vline(xintercept = c(1900, 2025), linetype = 2, colour = "grey40") + facet +
  labs(x = "Year", y = NULL) + theme_bw(base_size = 15) +
  theme(aspect.ratio = 1, strip.background = element_blank(), strip.text = element_text(hjust = 0), panel.grid.minor = element_blank())
p_all <- trajplot(ts, facet_grid(Trait ~ HT, scales = "free_y", labeller = label_parsed))
ggsave(file.path(out, "pred_growth_traj.pdf"), p_all, width = 10, height = 9); ggsave(file.path(out, "pred_growth_traj.png"), p_all, width = 10, height = 9, dpi = 160)
p_med <- trajplot(ts[HTarget == 0.1], facet_wrap(~Trait, scales = "free_y", labeller = label_parsed))
ggsave(file.path(out, "pred_growth_traj_median.pdf"), p_med, width = 12, height = 4.3); ggsave(file.path(out, "pred_growth_traj_median.png"), p_med, width = 12, height = 4.3, dpi = 180)
## ---- 2) per-run summaries: century 1900-2025 and projection 2025-2125 ----
smr <- function(d) { x <- d[G >= 10000 & G <= APOT]; at <- function(v, G0) d[G == G0, get(v)][1]
  rbindlist(lapply(c(LT = "LT", g = "g", hmax = "hmax"), function(v) data.table(trait = v,
    r_total = log(at(v, APOT) / at(v, 10000)), pct_total = 100 * (at(v, APOT) / at(v, 10000) - 1),
    pct_proj = 100 * (at(v, 10225) / at(v, APOT) - 1), level_2025 = at(v, APOT),
    h2_bar = mean(x[[paste0("Va_", v)]] / x[[paste0("Vp_", v)]]), CV_bar = 100 * mean(sqrt(x[[paste0("Va_", v)]]) / x[[v]])))) }
rs <- runs[, smr(dt[[1]]), by = .(file, Ve, sd, K, HTarget, Hx, w, f, SLT, Shmax, theta)]
g25 <- rs[trait == "g", .(file, g_2025 = level_2025)]; rs <- merge(rs, g25, by = "file"); rs[, wg := w * g_2025]
fwrite(rs, file.path(out, "pps_run_summaries.csv"))
q <- function(x) sprintf("%.2f [%.2f, %.2f]", median(x), quantile(x, .05), quantile(x, .95))
## headline: baseline, continued fishing (Hx1) vs fishing stopped (Hx0)
head_tab <- rs[eval(base), .(n = .N, change_1900_2025 = q(pct_total), change_2025_2125 = q(pct_proj), h2 = q(h2_bar), CV_A = q(CV_bar)), by = .(HTarget, Hx, trait)][order(trait, HTarget, -Hx)]
fwrite(head_tab, file.path(out, "pps_headline.csv")); print(head_tab)
## recovery: share of the 1900-2025 change reversed by 2125 after fishing stops
rec <- rs[eval(base) & Hx == 0, .(n = .N, recovered_pct_of_change = q(-100 * ((1 + pct_proj / 100) - 1) / ((1 + pct_total / 100) - 1 + 1e-12))), by = .(HTarget, trait)]
print(rec)
## ---- 3) sensitivity: violins at H = 0.10 (sd, K, Ve one at a time) ----
dv <- rs[Hx == 1 & HTarget == 0.1][, trait_lab := factor(trait, c("hmax", "g", "LT"), c("h[max]", "g", "L[T]"))]
panel <- function(xv, flt, xlab, yv, ylab, zero = FALSE) { d <- dv[eval(flt)]
  p <- ggplot(d, aes(factor(get(xv)), get(yv))) + geom_violin(fill = "grey85", colour = "grey40", trim = FALSE) +
    geom_boxplot(width = .18, outlier.size = .6, alpha = .6) + facet_grid(trait_lab ~ ., scales = "free_y", labeller = labeller(.rows = label_parsed)) +
    labs(x = xlab, y = ylab) + theme_bw() + theme(panel.grid.minor = element_blank(), strip.background = element_blank())
  if (zero) p <- p + geom_hline(yintercept = 0, linetype = 2, colour = "grey40"); p }
bsd <- quote(Ve == 0.05 & K == 10000); bK <- quote(Ve == 0.05 & sd == 0.05); bVe <- quote(K == 10000 & sd == 0.05)
yl <- "change 1900-2025 (%)"
pdf(file.path(out, "Scenarios_summary.pdf"), width = 12, height = 10)
print((panel("sd", bsd, "mutational SD (sd)", "pct_total", yl, TRUE) + ggtitle("a") | panel("K", bK, "carrying capacity (K)", "pct_total", yl, TRUE) + ggtitle("b") | panel("Ve", bVe, expression("environmental variance ("*V[e]*")"), "pct_total", yl, TRUE) + ggtitle("c")) /
  (panel("sd", bsd, "mutational SD (sd)", "CV_bar", expression(bar(CV)[A])) + ggtitle("d") | panel("K", bK, "carrying capacity (K)", "CV_bar", expression(bar(CV)[A])) + ggtitle("e") | panel("Ve", bVe, expression("environmental variance ("*V[e]*")"), "CV_bar", expression(bar(CV)[A])) + ggtitle("f")) /
  (panel("sd", bsd, "mutational SD (sd)", "h2_bar", expression(bar(h)^2)) + ggtitle("g") | panel("K", bK, "carrying capacity (K)", "h2_bar", expression(bar(h)^2)) + ggtitle("h") | panel("Ve", bVe, expression("environmental variance ("*V[e]*")"), "h2_bar", expression(bar(h)^2)) + ggtitle("i")))
dev.off()
sens <- rs[Hx == 1, .(n = .N, change_1900_2025 = q(pct_total), h2 = q(h2_bar), CV_A = q(CV_bar)), by = .(HTarget, trait, Ve, sd, K)][order(trait, HTarget, Ve, sd, K)]
fwrite(sens, file.path(out, "pps_sensitivity_by_scenario.csv"))
## ---- 4) eta^2 (all Hx1 runs, H as categorical control), as in R/PPS.R ----
df <- rs[Hx == 1][, `:=`(Ve_z = as.numeric(scale(Ve)), sd_z = as.numeric(scale(sd)), logK_z = as.numeric(scale(log(K))), H_fac = factor(HTarget))]
blk <- function(fit, model, keep, lab) { cf <- fixef(fit)[keep]; ci <- confint(fit, parm = keep, method = "Wald")
  et <- as.data.table(eta_squared(fit, partial = TRUE, ci = .95, alternative = "two"))[Parameter %in% keep]; r2s <- r2(fit)
  data.table(Model = model, term = keep, Term = lab[keep], beta = cf, ci_low = ci[, 1], ci_high = ci[, 2],
             eta2p = et$Eta2_partial[match(keep, et$Parameter)], eta2p_low = et$CI_low[match(keep, et$Parameter)], eta2p_high = et$CI_high[match(keep, et$Parameter)],
             R2_m = r2s$R2_marginal, R2_c = r2s$R2_conditional) }
lab <- c(Ve_z = "V[e]", sd_z = "sd", logK_z = "log(K)", h2_bar = "h^2", CV_bar = "CV[A]")
eta <- rbindlist(lapply(c("g", "LT", "hmax"), function(tr) { d <- df[trait == tr]
  rbind(blk(lmer(h2_bar ~ Ve_z + sd_z + logK_z + H_fac + (1 | theta), d), '"a)"~ h^2 %~% V[e] + sd + log(K)', c("Ve_z", "sd_z", "logK_z"), lab),
        blk(lmer(CV_bar ~ Ve_z + sd_z + logK_z + H_fac + (1 | theta), d), '"b)"~ CV[A] %~% V[e] + sd + log(K)', c("Ve_z", "sd_z", "logK_z"), lab),
        blk(lmer(r_total ~ h2_bar + CV_bar + logK_z + H_fac + (1 | theta), d), '"c)"~ r[T] %~% h^2 + CV[A] + log(K)', c("h2_bar", "CV_bar", "logK_z"), lab))[, trait := tr] }))
fwrite(eta, file.path(out, "eta_table.csv"))
eta[, Term := factor(Term, levels = c("V[e]", "sd", "log(K)", "h^2", "CV[A]"))]
pd <- position_dodge(width = .7)
p_eta <- ggplot(eta, aes(trait, eta2p, fill = Term)) + geom_col(position = pd, width = .6, colour = "grey25") +
  geom_errorbar(aes(ymin = eta2p_low, ymax = eta2p_high), position = pd, width = .15, colour = "grey20") +
  facet_grid(. ~ Model, labeller = labeller(.cols = label_parsed)) +
  scale_fill_manual(values = pal_jco()(10), labels = function(x) parse(text = x), name = "Predictor") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), name = expression(eta[partial]^2)) +
  scale_x_discrete(labels = function(x) parse(text = c(g = "g", hmax = "h[max]", LT = "L[T]")[x]), name = NULL) +
  theme_bw() + theme(aspect.ratio = 1, strip.text.x = element_text(hjust = 0), panel.grid.minor = element_blank(), strip.background = element_blank(), legend.position = "bottom")
ggsave(file.path(out, "eta_sum.pdf"), p_eta, width = 8, height = 3.5); ggsave(file.path(out, "eta_sum.png"), p_eta, width = 8, height = 3.5, dpi = 200)
## LaTeX table (replaces Tables/Table_Parametersizes_main.tex content; written here only)
tex <- eta[, sprintf("%s & %s & %s & %.3f [%.3f, %.3f] & %.3f [%.3f, %.3f] & %.3f & %.3f \\\\",
  c(h2_bar = "$\\overline{h^2}$", CV_bar = "$\\overline{CV_A}$", r_total = "$r_T$")[ifelse(grepl("^\"a", Model), "h2_bar", ifelse(grepl("^\"b", Model), "CV_bar", "r_total"))],
  c(Ve_z = "$V_e$", sd_z = "$sd$", logK_z = "$\\log K$", h2_bar = "$\\overline{h^2}$", CV_bar = "$\\overline{CV_A}$")[term],
  c(g = "$g$", LT = "$L_T$", hmax = "$h_{max}$")[trait], beta, ci_low, ci_high, eta2p, eta2p_low, eta2p_high, R2_m, R2_c)]
writeLines(c("\\begin{table}[h]", "\\caption{Standardized effects of environmental variance ($V_e$), mutational SD ($sd$) and carrying capacity ($\\log K$) on mean heritability, mean coefficient of additive genetic variation, and the log response $r_T$ (1900--2025), from the posterior predictive simulations (final posterior). Harvest level ($H = 0.02, 0.10, 0.16$) was included as a categorical control; its coefficients are not shown. Entries: standardized coefficient [95\\% CI], partial $\\eta^2$ [95\\% CI], marginal and conditional $R^2$.}",
  "\\label{Parametersizes_main}", "\\begin{tabular}{lllllll}", "\\toprule", "Response & Predictor & Trait & Std. $\\beta$ [95\\% CI] & Partial $\\eta^2$ [95\\% CI] & $R^2_m$ & $R^2_c$ \\\\", "\\midrule", tex, "\\bottomrule", "\\end{tabular}", "\\end{table}"),
  file.path(out, "Table_Parametersizes_main.tex"))
## ---- 5) post-stratified by w (and w*g), baseline scenarios ----
wb <- c(0.1, 0.5, 0.9, 1.2, 1.5)
ps <- rs[eval(base) & Hx == 1][, w_band := cut(w, wb, include.lowest = TRUE)]
ps_w <- ps[, .(n = .N, change_1900_2025 = q(pct_total), change_2025_2125 = q(pct_proj), median = median(pct_total)), by = .(HTarget, trait, w_band)][order(trait, HTarget, w_band)]
fwrite(ps_w, file.path(out, "pps_poststrat_w.csv"))
ps[, wg_tercile := cut(wg, quantile(wg, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, dig.lab = 2), by = HTarget]
ps_wg <- ps[, .(n = .N, change_1900_2025 = q(pct_total), median = median(pct_total)), by = .(HTarget, trait, wg_tercile)][order(trait, HTarget, wg_tercile)]
fwrite(ps_wg, file.path(out, "pps_poststrat_wg.csv"))
print(ps_w[trait == "LT"]); print(ps_wg[trait == "LT"])
p_w <- ggplot(ps[trait == "LT"], aes(w, pct_total, colour = factor(HTarget))) + geom_hline(yintercept = 0, linetype = 3) +
  geom_point(alpha = .5, size = 1.2) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  scale_colour_viridis_d("scenario H", end = .85) + theme_bw(base_size = 13) + theme(panel.grid.minor = element_blank()) +
  labs(x = "cost of reproduction w", y = "change in length at maturation L_T, 1900-2025 (%)")
ggsave(file.path(out, "pps_dLT_vs_w.pdf"), p_w, width = 7, height = 5); ggsave(file.path(out, "pps_dLT_vs_w.png"), p_w, width = 7, height = 5, dpi = 180)
