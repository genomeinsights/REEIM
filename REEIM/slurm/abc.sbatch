#!/bin/bash
#SBATCH -J reeim-abc
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p main
#SBATCH -c 128
#SBATCH -t 24:00:00
#SBATCH -o slurm-%x-%j.out
#
# ABC-SMC for one study system. Submit from the repository root:
#   sbatch -A <project> [-p <partition>] [-c <cores>] [-t <time>] REEIM/slurm/abc.sbatch perch
# The system name picks the input folder REEIM/input/<system>; results go to REEIM/out_<system>/.
# If the job hits its time limit, submit the same command again: it continues where it stopped.
set -eo pipefail
SYSTEM=${1:-perch}
cd "${SLURM_SUBMIT_DIR:-.}"
source REEIM/slurm/env.sh
INPUT=./REEIM/input/$SYSTEM
OUT=./REEIM/out_$SYSTEM/
[ -d "$INPUT" ] || { echo "no input folder $INPUT"; exit 1; }
mkdir -p "$OUT"
# last cycle compared with the data = burn-in + fishing years - 1 (as in run_REEIM_ABC.R)
TDELTA=$(awk -F, '$1=="Tdelta"{print $2}' "$INPUT/parameters.csv")
LAST=$((TDELTA + $(tail -n +2 "$INPUT/fishing_intensity.csv" | grep -c .) - 1))
# remove runs cut off by an earlier time limit, in the generation that was in progress
CUR=$(ls -d "$OUT"GEN_* 2>/dev/null | sort -V | tail -n 1 || true)
[ -n "$CUR" ] && bash REEIM/slurm/clean_partial.sh "$CUR" "$LAST"
echo "system $SYSTEM, $SLURM_CPUS_PER_TASK cores, slim $(slim -v | head -n 1), commit $(git rev-parse --short HEAD 2>/dev/null)"
Rscript REEIM/run_REEIM_ABC.R "$PWD" "${SLURM_CPUS_PER_TASK:-4}" "$INPUT" "$OUT"
