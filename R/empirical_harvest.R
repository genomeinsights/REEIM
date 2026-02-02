HRate <- fread("./empirical_data/harvestrate.csv")

# Accepted particles from final generation
base_folder <- "./REEIM/out_REEIM/"
APOT        <- 10125     
REEIM_res   <- readRDS(file.path(base_folder, "REEIM_res.rds"))
particles   <- REEIM_res$particles
post        <- particles[[length(particles)]][keep == TRUE]


qq_data <- data.table(do.call(cbind,qqplot(post[which(keep)]$H*100,HRate[,harvestperc])))
setnames(qq_data,c("x","y"), c("Posterior dist. for H (%)","Empirical estim. (%)"))
plot(qq_data,ylim=c(0,25),xlim=c(0,25))
abline(lm(`Empirical estim. (%)`~`Posterior dist. for H (%)`+0,data=qq_data[`Posterior dist. for H (%)`<15]),lty=2)
abline(0,1)
text(23,3,"slope=0.27")



qq_data[`Posterior dist. for H (%)`<15,cor(`Posterior dist. for H (%)`,`Empirical estim. (%)`,method = "spearman")^2]

qq_data[, `Posterior rank` := rank(`Posterior dist. for H (%)`) / .N]
qq_data[, `Empirical rank` := rank(`Empirical estim. (%)`) / .N]

plot(
  qq_data$`Posterior rank`,
  qq_data$`Empirical rank`,
  xlab = "Posterior rank (H)",
  ylab = "Empirical rank",
  pch = 16
)
abline(0,1,lty=2)
