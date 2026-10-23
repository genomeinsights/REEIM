############################################################
##  In-silico recovery test, step 1: extract per-draw stats
##
##  Reads the generation-1 (prior) ABC draws written by
##  REEIM/run_REEIM_ABC.R and stores, for every completed draw,
##  everything needed to use it as a pseudo-observed dataset:
##    - the parameters (parsed from the file name, as in the driver)
##    - the size-density curve at APOT (same KDE as the ABC distance)
##    - the mean size-at-age curve at APOT
##    - L_T at APOT (the third ABC summary)
##    - the century trait responses r_T (t0 = 10000, t1 = APOT),
##      computed as in R/trati_cor_with_H.R
##
##  Run from the repository root:
##    Rscript recovery_test/01_extract_gen1_stats.R [gen_folder] [cores]
##  Output: recovery_test/output/gen1_stats.rds
############################################################

library(data.table)
library(parallel)

args       <- commandArgs(trailingOnly = TRUE)
gen_folder <- if (length(args) > 0) args[1] else "./REEIM/out_REEIM/GEN_1"
cores      <- if (length(args) > 1) as.integer(args[2]) else max(1, detectCores() - 1)
res_file   <- "./REEIM/out_REEIM/REEIM_res.rds"
out_dir    <- "./recovery_test/output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

source("./recovery_test/R/functions.R")

## Which draws to use: the 2,500 kept in REEIM_res if it exists, else all
files <- list.files(gen_folder, full.names = TRUE, pattern = "^demographicData_")
if (file.exists(res_file)) {
  kept  <- readRDS(res_file)$particles[[1]]$file
  files <- files[gsub("demographicData_", "", basename(files)) %in% kept]
  message(length(files), " generation-1 draws matched REEIM_res.rds")
} else {
  message("REEIM_res.rds not found, using all ", length(files), " draws in ", gen_folder)
}

res  <- mclapply(files, extract_one, mc.cores = cores)
gen1 <- bind_draws(res)
message(nrow(gen1$theta), " draws extracted")
saveRDS(gen1, file.path(out_dir, "gen1_stats.rds"))
