#!/bin/bash
# Usage: REEIM/slurm/clean_partial.sh <output folder> <last cycle>
# Deletes runs that were cut off by a job time limit: a traitStats file that
# doesn't reach <last cycle>, together with its demographicData file.
# Only call this when no SLiM job is writing to the folder.
dir=$1; last=$2; n=0
for t in "$dir"/traitStats_*.txt; do
  [ -e "$t" ] || continue
  g=$(tail -n 1 "$t" | awk '{print $1}')
  if ! [[ "$g" =~ ^[0-9]+$ ]] || [ "$g" -lt "$last" ]; then
    rm -f "$t" "${t/traitStats_/demographicData_}"; n=$((n+1))
  fi
done
for d in "$dir"/demographicData_*.txt; do
  [ -e "$d" ] || continue
  [ -e "${d/demographicData_/traitStats_}" ] || { rm -f "$d"; n=$((n+1)); }
done
echo "clean_partial: removed $n unfinished runs in $dir"
