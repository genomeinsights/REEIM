library(data.table)
library(readxl)

## Size data for the size-distribution comparison in the ABC:
##   "nors"       NORS export of 2025-06-06 (nors250606.csv): perch, 2010-2024, 62 lakes, lengths >= 6 cm
##   "monitoring" Monitoring_data.csv: perch, 2010-2020, 79 lakes (used for the published fit)
size_source <- "nors"

age_structure_obs <- as.data.table(read_excel("./empirical_data/ÅlderSjöprovfiskenAbborre.xlsx"))[,.(age=Ålder,size=Längd)]
size_structure_obs <- if (size_source == "nors") {
  fread("./empirical_data/Voulnerability/nors250606.csv")[, .(Lake_ID = vatten_id, Year, L = cmklass)]
} else {
  fread("./empirical_data/Voulnerability/Monitoring_data.csv")[SpeciesName=="perch",.(Lake_ID,  Year, L=cmklass)]
}

## Replace only these two tables; other objects already in the file (e.g. growth_curve_obs) are kept
out_file <- "./empirical_data/observed_data.RData"
obs <- new.env()
if (file.exists(out_file)) load(out_file, envir = obs)
assign("age_structure_obs", age_structure_obs, envir = obs)
assign("size_structure_obs", size_structure_obs, envir = obs)
save(list = ls(obs), envir = obs, file = out_file)
