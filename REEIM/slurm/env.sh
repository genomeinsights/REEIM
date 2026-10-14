# Sourced by the job scripts: make slim, Rscript and GNU parallel available.
# Edit this file for your cluster if the conda setup differs.
if ! command -v conda >/dev/null 2>&1; then module load conda 2>/dev/null || true; fi
eval "$(conda shell.bash hook)"
conda activate reeim
