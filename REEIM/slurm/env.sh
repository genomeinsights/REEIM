# Sourced by the job scripts: make slim, Rscript and GNU parallel available.
# Edit this file for your cluster if the conda setup differs.
if ! command -v conda >/dev/null 2>&1; then module load PDC/26.03 miniconda3/26.3.2-2-cpeGNU-26.03 2>/dev/null || module load conda 2>/dev/null || true; fi   # Dardel, else generic
eval "$(conda shell.bash hook)"
conda activate reeim
