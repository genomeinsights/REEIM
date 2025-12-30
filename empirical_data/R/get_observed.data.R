library(data.table)
library(ggplot2)
library(readxl)


age_structure_obs <- as.data.table(read_excel("./empirical_data/ÅlderSjöprovfiskenAbborre.xlsx"))[,.(age=Ålder,size=Längd)]
size_structure_obs <- fread("./empirical_data/voulnerability/Monitoring_data.csv")[SpeciesName=="perch",.(Lake_ID,  Year, L=cmklass)]
save(age_structure_obs,size_structure_obs,file="./empirical_data/observed_data.RData")

