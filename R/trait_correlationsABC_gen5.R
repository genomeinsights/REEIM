library(data.table)
library(egg)
library(ggplot2)
library(viridis)

# Accepted particles from final generation
base_folder <- "./REEIM/out_REEIM/"
APOT        <- 10125   
REEIM_res <- readRDS(file.path(base_folder, "REEIM_res.rds"))
particles   <- REEIM_res$particles
post <- particles[[length(particles)]][keep == TRUE]

# Keep only key parameters
pp <- post[, .(H, h, LT, g, weights)]

#pp[,cor.test(H,h)]

p_H_h <- ggplot(pp, aes(H, h)) +
  geom_point(
    aes(alpha = weights),
    color = "grey30",
    size = 1
  ) +
  geom_smooth(
    aes(weight = weights),
    method = "lm",
    se = FALSE,
    color = "firebrick",
    linewidth = 1
  ) +
  scale_alpha(range = c(0.2, 1), guide = "none") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 0.5)+
  labs(
    x = expression(H),
    y = expression(h),
    title = "a) Harvest intensity and compensatory growth"
  )

p_H_LT <- ggplot(pp, aes(H, LT)) +
  geom_point(
    aes(alpha = weights),
    color = "grey30",
    size = 1
  ) +
  geom_smooth(
    aes(weight = weights),
    method = "lm",
    se = FALSE,
    color = "firebrick",
    linewidth = 1
  ) +
  scale_alpha(range = c(0.2, 1), guide = "none") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 0.5)+
  labs(
    x = expression(H),
    y = expression(L[T]),
    title = "b) Harvest intensity and maturation size"
  )
p_H_g <- ggplot(pp, aes(H, g)) +
  geom_point(
    aes(alpha = weights),
    color = "grey30",
    size = 1
  ) +
  geom_smooth(
    aes(weight = weights),
    method = "lm",
    se = FALSE,
    color = "firebrick",
    linewidth = 1
  ) +
  theme(aspect.ratio = 0.5)+
  scale_alpha(range = c(0.2, 1), guide = "none") +
  theme_bw(base_size = 14) +
  labs(
    x = expression(H),
    y = expression(g),
    title = "c) Harvest intensity and parental investment"
  )


fig_corr <- ggarrange(
  p_H_h,
  p_H_LT,
  p_H_g,
  ncol = 3
)

fig_corr

ggsave(
  "./Figures/par_cor.pdf",
  fig_corr,
  width = 18,
  height = 7
)