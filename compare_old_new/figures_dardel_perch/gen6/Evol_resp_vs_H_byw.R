library(data.table); library(ggplot2); t0 <- 10000; t1 <- 10125
a <- commandArgs(TRUE); base_folder <- a[1]; GEN <- as.integer(a[2]); out <- a[3]; bands <- as.numeric(a[4:5])
APOT        <- 10125   
REEIM_res <- readRDS(file.path(base_folder, "REEIM_res.rds"))
acc   <- REEIM_res$particles[[GEN]][which(keep)]

run_root <- file.path(base_folder, paste0("GEN_", GEN))
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
  s[, `:=`(H = acc$H[i], w = acc$w[i])]
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
dat_corr <- run_summaries_base[!is.na(H) & !is.na(pct_total)]
## w bands: given cut-offs (low < bands[1], high > bands[2]); NA -> weighted terciles of w
if (anyNA(bands)) { o <- acc[order(w)]; cw <- cumsum(o$weights) / sum(o$weights); bands <- c(o$w[which(cw >= 1/3)[1]], o$w[which(cw >= 2/3)[1]]) }
dat_corr[, band := fifelse(w < bands[1], "low", fifelse(w > bands[2], "high", "mid"))]
band_lab <- c(low = sprintf("w < %.2f", bands[1]), high = sprintf("w > %.2f", bands[2]))

fit_tbl <- function(d) { fitw <- lm(r_total ~ H, weights = weight, data = d); ci <- confint(fitw)["H", ]
  wn <- d$weight / sum(d$weight); Hp <- d$H - sum(wn * d$H); Rp <- d$r_total - sum(wn * d$r_total)
  data.table(n = nrow(d), slope_H = coef(fitw)[["H"]], slope_ci_lo = ci[1], slope_ci_hi = ci[2],
             pearson_r = sum(wn * Hp * Rp) / sqrt(sum(wn * Hp^2) * sum(wn * Rp^2)),
             pct_per_0.1H = 100 * (exp(coef(fitw)[["H"]] / 10) - 1)) }
w_tbl   <- dat_corr[, fit_tbl(.SD), by = trait_plot]
band_tbl <- dat_corr[band != "mid", fit_tbl(.SD), by = .(trait_plot, band)][, band := band_lab[band]][order(trait_plot, band)]
cat(sprintf("w band cut-offs: %.3f / %.3f\n", bands[1], bands[2])); print(w_tbl); print(band_tbl)

line_dt <- dat_corr[band != "mid", { f <- lm(pct_total ~ H, weights = weight); h <- range(H); .(H = h, pct_total = predict(f, data.frame(H = h))) }, by = .(trait_plot, band)][, band := factor(band_lab[band], levels = band_lab)]
corr_annot_w <- w_tbl[, .(trait_plot, label = sprintf("italic(r)==%.2f~~beta[H]==%.3f", pearson_r, slope_H))]

mk <- function(lines) {
  p <- ggplot(dat_corr[order(w)], aes(x = H, y = pct_total)) +
    geom_point(aes(colour = w), alpha = 0.6, size = 1.8) +
    geom_smooth(aes(group = trait_plot), method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.8, col = "black") +
    scale_colour_viridis_c("cost of\nreproduction w") +
    facet_grid(. ~ trait_plot, labeller = label_parsed) +
    scale_x_continuous("harvest intensity H (fraction of biomass removed per year)", limits = c(0, NA)) +
    scale_y_continuous(expression(r[T]~", %"), labels = function(v) sprintf("%.1f", v)) +
    geom_text(data = corr_annot_w, aes(x = 0.003, y = Inf, label = label), parse = TRUE, hjust = 0, vjust = 1.5, size = 3.5) +
    theme_bw(base_size = 11) + theme(aspect.ratio = 1, panel.grid.minor = element_blank(), strip.background = element_blank())
  if (lines) p <- p + geom_line(data = line_dt, aes(linetype = band), colour = "grey25", linewidth = 0.45) +
    scale_linetype_manual(NULL, values = c(2, 3))
  p }
for (L in c(TRUE, FALSE)) { s <- if (L) "" else "_nolines"
  ggsave(paste0(out, s, ".pdf"), mk(L), height = 4, width = 13); ggsave(paste0(out, s, ".png"), mk(L), height = 4, width = 13, dpi = 200) }
fwrite(band_tbl, paste0(out, "_band_slopes.csv"))
