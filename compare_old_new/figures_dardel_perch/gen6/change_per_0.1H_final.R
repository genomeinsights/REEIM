## Final per-0.1 H responses, same estimator as the Fig. 4 beta labels:
## weighted lm of log change (log X_2025 - log X_1900) on H; per 0.1 H = exp(0.1 beta) - 1 (CI converted the same way).
## Absolute g: weighted lm of (g_2025 - g_1900) on H, times 0.1.
suppressMessages(library(data.table)); a <- commandArgs(TRUE); r <- readRDS(a[1])
k <- r$particles[[length(r$particles)]][keep == TRUE]
tr <- fread(a[2], header = FALSE, col.names = c("tf", "LT0", "LT1", "g0", "g1", "h0", "h1"))[, file := sub("^traitStats_", "", tf)]
k <- merge(k, tr, by = "file")
f <- function(y, q, rel) { m <- lm(y ~ H, weights = weights, data = k); b <- coef(m)[["H"]]; ci <- confint(m)["H", ]
  tr <- if (rel) function(x) 100 * (exp(0.1 * x) - 1) else function(x) 0.1 * x
  data.table(quantity = q, beta_H = b, per_0.1_H = tr(b), ci_lo = tr(ci[1]), ci_hi = tr(ci[2]),
             unit = if (rel) "% change" else "absolute change in g", n = nrow(k)) }
out <- rbind(f(log(k$LT1 / k$LT0), "L_T", TRUE), f(log(k$h1 / k$h0), "h_max", TRUE), f(log(k$g1 / k$g0), "g", TRUE), f(k$g1 - k$g0, "g_abs", FALSE))
print(out, digits = 4); fwrite(out, a[3])
