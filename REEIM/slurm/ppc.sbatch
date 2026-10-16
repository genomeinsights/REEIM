#!/bin/bash
#SBATCH -J reeim-ppc
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p main
#SBATCH -c 256
#SBATCH -t 06:00:00
#SBATCH -o slurm-%x-%j.out
#
# Posterior predictive check of evolutionary change for a plaice ABC, after that ABC has finished.
# Submit from the repository root:
#   sbatch -A <project> [-p <partition>] [-c <cores>] REEIM/slurm/ppc.sbatch plaice [n_draws=250]
# REEIM_ABC_OUT / REEIM_PPC_OUT override the input and output folders (test runs); REEIM_SEED (default 1) the draws.
# Reads REEIM/out_<system>/REEIM_res.rds and writes REEIM/ppc_<system>/:
#   ppc_genetic_summary.csv  (genetic L_T change over the accepted particles; no new runs)
#   ppc_survey_summary.csv   (genetic L_T, survey ogive L50 and survey PMRN Lp50 age 4 over n_draws new runs
#                             of REEIM_vesurvey.slim) and fig/ppc_trajectories.png.
#   ppc_sensitivity_summary.csv (fishing effect on L_T per 0.1 of H, from n_draws paired H = 0 control runs,
#                             so 2 x n_draws runs in all: 500 runs, about 1 node-hour at 256 runs at once).
# Resubmitting skips the draws whose runs have finished.
set -eo pipefail
SYSTEM=${1:-plaice}; N=${2:-250}
cd "${SLURM_SUBMIT_DIR:-.}"
source REEIM/slurm/env.sh
ABC=${REEIM_ABC_OUT:-./REEIM/out_$SYSTEM/}; OUT=${REEIM_PPC_OUT:-./REEIM/ppc_$SYSTEM/}; INPUT=./REEIM/input/$SYSTEM   # overrides for test runs
SEED=${REEIM_SEED:-1}
[ -f "$ABC/REEIM_res.rds" ] || { echo "no $ABC/REEIM_res.rds"; exit 1; }
echo "system $SYSTEM, $N draws, ${SLURM_CPUS_PER_TASK:-4} cores, slim $(slim -v | head -n 1), commit $(git rev-parse --short HEAD 2>/dev/null)"
Rscript REEIM/plaice_validation/ppc_plaice.R genetic "$ABC" "$OUT"
[ -f "$OUT/cmds.txt" ] || Rscript REEIM/plaice_validation/ppc_plaice.R prepare "$ABC" "$OUT" "$INPUT" "$N" "$SEED"
# Run manifest (as abc.sbatch): code commit and input checksums once per submission; cmds.txt holds each run's
# seed and full command, joblog.txt each run's exit status.
{ echo "# commit $(git rev-parse HEAD 2>/dev/null)"; echo "# date $(date -u +%FT%TZ)"; echo "# slim $(slim -v | head -n 1)"
  md5sum "$INPUT"/*.csv REEIM/plaice_validation/test_runs/REEIM_vesurvey.slim REEIM/plaice_validation/ppc_plaice.R "$ABC/REEIM_res.rds" 2>/dev/null | sed 's/^/# md5 /'
} > "$OUT/manifest_header_$(date +%Y%m%d_%H%M%S).txt"

# run only the draws without a finished traitStats file (reaching year 2020 = cycle Tdelta + 170)
TDELTA=$(awk -F, '$1=="Tdelta"{print $2}' "$INPUT/parameters.csv"); LAST=$((TDELTA + 170))
: > "$OUT/todo.txt"
while IFS= read -r c; do
  r=$(echo "$c" | sed -E 's/.* -d rep=([0-9]+) .*/\1/')
  f=$(ls "$OUT"/out/traitStats_*_rep"$r".txt 2>/dev/null | head -n 1 || true)
  g=$( [ -n "$f" ] && tail -n 1 "$f" | cut -d' ' -f1 | grep -E '^[0-9]+$' || echo 0 )
  [ "$g" -ge "$LAST" ] || echo "$c" >> "$OUT/todo.txt"
done < "$OUT/cmds.txt"
echo "$(grep -c . "$OUT/todo.txt" || true) runs to do"
parallel --joblog "+$OUT/joblog.txt" -j "${SLURM_CPUS_PER_TASK:-4}" < "$OUT/todo.txt" || true
Rscript REEIM/plaice_validation/ppc_plaice.R survey "$OUT"
