#!/bin/bash
#SBATCH -J reeim-pps
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p main
#SBATCH -c 128
#SBATCH -t 24:00:00
#SBATCH -o slurm-%x-%j.out
#
# Posterior predictive simulations (perch), after the perch ABC has finished. From the repository root:
#   sbatch -A <project> [-p <partition>] REEIM/slurm/pps.sbatch
# Reads REEIM/out_perch/REEIM_res.rds and writes REEIM/scenarios_perch/.
# Resubmitting skips scenarios whose folder exists, so if a job is cut off,
# delete the last (incomplete) scenario folder before resubmitting.
set -eo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source REEIM/slurm/env.sh
export REEIM_CORES=${SLURM_CPUS_PER_TASK:-4}
export REEIM_ABC_OUT=./REEIM/out_perch/
export REEIM_PPS_OUT=./REEIM/scenarios_perch/
echo "$REEIM_CORES cores, slim $(slim -v | head -n 1), commit $(git rev-parse --short HEAD 2>/dev/null)"
Rscript REEIM/RUN_PPS.R
