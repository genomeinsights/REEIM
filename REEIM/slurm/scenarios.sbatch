#!/bin/bash
#SBATCH -J reeim-scen
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p main
#SBATCH -c 128
#SBATCH -t 12:00:00
#SBATCH -o slurm-%x-%j.out
#
# Evolutionary-change scenarios for a plaice ABC (after that ABC has finished): posterior draws x genetics x H.
# Submit from the repository root:
#   sbatch -A <project> [-c <cores>] REEIM/slurm/scenarios.sbatch plaice [n_draws=50] [H list=0,0.04,0.078,0.12,0.2]
# Genetics: REEIM/plaice_validation/scenario_genetics.csv (Ve fixed; mutation rate u sets Va, since Va scales with u
# and hardly with sd in REEIM; realised h2 ~ 0.1-0.6). Writes REEIM/scen_<system>/: scenario_summary.csv,
# scenario_runs.csv, fig/scenario_response.png and fig/scenario_fit.png. 50 draws x 5 x 5 = 1250 runs, ~2-3 h on one
# node (high-u runs are slower). Resubmitting skips finished runs.
set -eo pipefail
SYSTEM=${1:-plaice}; N=${2:-50}; HS=${3:-0,0.04,0.078,0.12,0.2}
cd "${SLURM_SUBMIT_DIR:-.}"
source REEIM/slurm/env.sh
ABC=./REEIM/out_$SYSTEM/; OUT=./REEIM/scen_$SYSTEM/; INPUT=./REEIM/input/$SYSTEM
[ -f "$ABC/REEIM_res.rds" ] || { echo "no $ABC/REEIM_res.rds"; exit 1; }
echo "system $SYSTEM, $N draws, H $HS, ${SLURM_CPUS_PER_TASK:-4} cores, slim $(slim -v | head -n 1), commit $(git rev-parse --short HEAD 2>/dev/null)"
[ -f "$OUT/cmds.txt" ] || Rscript REEIM/plaice_validation/ppc_plaice.R scenarios "$ABC" "$OUT" "$INPUT" REEIM/plaice_validation/scenario_genetics.csv "$N" "$HS"
TDELTA=$(awk -F, '$1=="Tdelta"{print $2}' "$INPUT/parameters.csv"); LAST=$((TDELTA + 170))
: > "$OUT/todo.txt"
while IFS= read -r c; do
  r=$(echo "$c" | sed -E 's/.* -d rep=([0-9]+) .*/\1/')
  f=$(ls "$OUT"/out/traitStats_*_rep"$r".txt 2>/dev/null | head -n 1 || true)
  g=$( [ -n "$f" ] && tail -n 1 "$f" | cut -d' ' -f1 | grep -E '^[0-9]+$' || echo 0 )
  [ "$g" -ge "$LAST" ] || echo "$c" >> "$OUT/todo.txt"
done < "$OUT/cmds.txt"
echo "$(grep -c . "$OUT/todo.txt" || true) runs to do"
parallel -j "${SLURM_CPUS_PER_TASK:-4}" < "$OUT/todo.txt" || true
Rscript REEIM/plaice_validation/ppc_plaice.R scenario_summary "$OUT"
