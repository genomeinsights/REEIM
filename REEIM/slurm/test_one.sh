#!/bin/bash
# Quick check of the setup: one SLiM run with typical perch values (a few minutes).
# Run from the repository root, on a compute node if login nodes don't allow it:
#   bash REEIM/slurm/test_one.sh
set -eo pipefail
source REEIM/slurm/env.sh
mkdir -p REEIM/slurm_test
time slim -s 1 -d "folder=\"$PWD/REEIM/slurm_test\"" -d "inputDir=\"$PWD/REEIM/input/perch\"" \
  -d H=0.05 -d f=0.76 -d w=0.41 -d Shmax=88.3 -d SLT=191 -d projection=0 REEIM/REEIM.slim > REEIM/slurm_test/log.txt 2>&1 || true
ls REEIM/slurm_test
tail -n 1 REEIM/slurm_test/traitStats_*.txt
