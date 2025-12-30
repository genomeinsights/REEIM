# ============================================================
# Fig. S*: Historical fishing intensity and size selectivity
# ============================================================

library(data.table)
library(ggplot2)
library(mgcv)
library(egg)
library(scales)

# -------------------------------
# Panel a: Historical fishing intensity
# -------------------------------

sportfiske_data <- data.table(
  Year = c(1939,1947,1950,1963,1975,1980,1990,2000,2006,2008,
           2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023),
  Fishermen = c(333334,666667,1333334,2000000,2000000,2000000,
                2200000,2300000,1000000,1000000,1700000,1600000,
                1600000,1400000,1400000,1300000,1600000,1700000,
                1500000,1200000,1300000)
)

# GAM smoother (used only to define relative trajectory)
mod <- gam(Fishermen ~ s(Year, bs = "tp", k = 12),
           family = nb(), data = sportfiske_data)

pred_years <- data.table(Year = 1900:2025)
pred_years[, Fishermen := predict(mod, pred_years, type = "response")]

pred_years$Fishermen <- as.numeric(formatC(signif(pred_years$Fishermen,digits=4), digits=4,format="fg"))

# Truncate extrapolation beyond last observation
pred_years[1900:2025>2010,Fishermen := pred_years[1900:2025==2010,Fishermen]]

pred_years[,Fishermen/max(Fishermen)]
p_hist <- ggplot() +
  geom_point(
    data = sportfiske_data,
    aes(Year, Fishermen/1e6),
    color = "grey30",
    size = 2.2
  ) +
  geom_line(
    data = pred_years,
    aes(Year, Fishermen/1e6),
    color = "firebrick",
    linewidth = 0.9
  ) +
  scale_y_continuous(labels = comma) +
  theme_bw(base_size = 13) +
  theme(
    aspect.ratio = 1,
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  ) +
  labs(
    x = "Year",
    y = "Number of recreational anglers (Million)",
    title = "b) Historical fishing intensity"
  )

# -------------------------------
# Panel b: Size-dependent vulnerability
# -------------------------------

vulnerability_clean <- fread(
  "./empirical_data/Voulnerability/vulnerability_clean.csv"
)

p_vuln <- ggplot(
  vulnerability_clean[size < 36],
  aes(size, scaled_vulnerability)
) +
  geom_line(color = "firebrick", linewidth = 0.9) +
  geom_point(color = "grey30", size = 1.6) +
  theme_bw(base_size = 13) +
  theme(
    aspect.ratio = 1,
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  ) +
  labs(
    x = "Body length (cm)",
    y = "Relative fishing vulnerability",
    title = "a) Size-selective vulnerability"
  )

# -------------------------------
# Combine panels
# -------------------------------

fig_hist_fish <- ggarrange(
  p_vuln,
  p_hist,
  ncol = 2,
  widths = c(1,1)
)

ggsave(
  "./Figures/Hist_fishing.pdf",
  fig_hist_fish,
  width  = 10,
  height = 5
)

SSF <- vulnerability_clean[size < 36,round(scaled_vulnerability,3)]
fix(SSF) ## copy this to REEIM.slim

FishingIntensity <- pred_years[,round(Fishermen/max(Fishermen),3)]
fix(FishingIntensity) ## copy this to REEIM.slim



