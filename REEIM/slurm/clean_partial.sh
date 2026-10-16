#!/bin/bash
# Usage: REEIM/slurm/clean_partial.sh <output folder> <last cycle>
# Deletes runs that were cut off by a job time limit: a run is kept only if its
# traitStats file has a row for <last cycle> and its demographicData file has the
# size counts for that cycle (written later in the same cycle). run_REEIM_ABC.R
# makes the same check itself; this script is a safety net before a restart.
# Only call this when no SLiM job is writing to the folder.
dir=$1; last=$2; n=0
for t in "$dir"/traitStats_*.txt; do
  [ -e "$t" ] || continue
  d=${t/traitStats_/demographicData_}
  if ! grep -q "^$last " "$t" || ! [ -e "$d" ] || ! grep -q "^countSize $last " "$d"; then
    rm -f "$t" "$d"; n=$((n+1))
  fi
done
for d in "$dir"/demographicData_*.txt; do
  [ -e "$d" ] || continue
  [ -e "${d/demographicData_/traitStats_}" ] || { rm -f "$d"; n=$((n+1)); }
done
echo "clean_partial: removed $n unfinished runs in $dir"
