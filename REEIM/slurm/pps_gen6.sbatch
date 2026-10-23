#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p main
#SBATCH -c 256
#SBATCH -t 14:00:00
#SBATCH -o slurm-%x-%j.out
# Perch PPS from the final posterior (GEN_6 of out_perch_final). Submit from the repository root with
#   sbatch -A <project> -J reeim-pps6-s<i> --export=ALL,REEIM_PPS_SHARD=<i>/<n> REEIM/slurm/pps_gen6.sbatch
set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source REEIM/slurm/env.sh
export REEIM_CORES=${SLURM_CPUS_PER_TASK:-4}
export REEIM_ABC_OUT=./REEIM/out_perch_final/
export REEIM_PPS_OUT=./REEIM/pps_perch_final_GEN6/
export REEIM_PPS_H=0.02,0.10,0.16
echo "$REEIM_CORES cores, shard ${REEIM_PPS_SHARD:-1/1}, slim $(slim -v | head -n 1), commit $(git rev-parse --short HEAD 2>/dev/null) (RUN_PPS.R modified, uncommitted)"
Rscript REEIM/RUN_PPS.R
