# Running the REEIM reruns on a SLURM cluster (Dardel, PDC)

These notes are for running the REEIM ABC-SMC reruns on a cluster that uses the SLURM scheduler.
Everything runs from one clone of the repository.
Each job uses the cores of one node, and GNU parallel spreads the SLiM runs over them.

## What is run, and in what order

| Step | Job script | Needs | Approx. size |
|---|---|---|---|
| 1. Perch ABC-SMC | `abc.sbatch perch` | nothing | 5 generations × 2,500 SLiM runs ≈ 1,100–1,700 core-hours |
| 2. Perch posterior predictive simulations | `pps.sbatch` | step 1 finished | 24 scenarios × 240 runs ≈ 500–800 core-hours |
| 3. Plaice ABC-SMC (two jobs) | `abc.sbatch plaice` (H prior 0.06–0.11, ICES catch/TSB) and `abc.sbatch plaice_Hfree` (H prior 0–0.2) | go-ahead from Petri (inputs still being tested) | similar to step 1, each |
| 3b. Plaice posterior predictive check | `ppc.sbatch plaice` and `ppc.sbatch plaice_Hfree` | the matching step-3 job finished | 250 runs each (~30 min on one node) |
| 4. Recovery test | (added later) | step 1; the scripts are rebuilt after it | about 30,000 runs |

The sizes assume about 5.5 min of core time per SLiM run, measured on a Mac mini running 13 at once, and up to 1.5× that on cluster cores.
Check with `test_one.sh` first.

**Dardel specifics.** A thin node has 128 physical cores (256 logical CPUs) and 256 GB of memory.
The `main` partition gives whole nodes, is charged for all 128 cores, and allows at most 24 hours per job.
The job scripts therefore default to `-p main -c 128 -t 24:00:00`, and one ABC generation (2,500 runs) takes about 2–3 hours on one node.
For smaller jobs, `-p shared` charges per core and allows up to 7 days.
Steps 1 and 3 (and the two plaice jobs) can run at the same time, because each system writes to its own folder.

## 1. Get the code

```
git clone https://github.com/genomeinsights/REEIM.git
cd REEIM
git checkout <commit given by Petri>      # branch fix/abc-rerun
```

Run every command below from the repository root.

## 2. Install the software (once)

The setup needs SLiM 5.0, R (with data.table, ks, mvtnorm and ggplot2) and GNU parallel.
All of them come from conda-forge:

```
module load conda          # or however conda is provided on the cluster
conda env create -f REEIM/slurm/environment.yml
```

`REEIM/slurm/env.sh` activates this environment in every job.
Edit that one file if your cluster provides conda differently.
Use SLiM 5.0: with SLiM 4.x each run exits with an error message after its last year.
The results are still written, but the log fills with errors.

## 3. Check one run

```
bash REEIM/slurm/test_one.sh        # on a compute node, e.g. inside an interactive job
```

The run takes a few minutes. Expect the `real` time to be roughly 3–6 minutes.
The last line printed should start with 10126, the final year.
Send Petri the time, because it sets the time limit to request.

## 4. Submit

```
sbatch -A <project> REEIM/slurm/abc.sbatch perch
```

- The defaults are one whole `main` node (128 cores) for 24 hours. Override them on the command line if needed, e.g. `-p shared -c 32` for a small test.
- To queue a continuation automatically, submit a second job that starts when the first ends: `sbatch -A <project> --dependency=afterany:<jobid> REEIM/slurm/abc.sbatch perch`.
- The job writes `REEIM/out_perch/GEN_1 … GEN_5` and, after each generation, `REEIM/out_perch/REEIM_res.rds`.
- Every SLiM run gets an explicit seed. `REEIM/out_perch/manifest.tsv` lists each run's generation, seed and full command, `manifest_header_*.txt` records the commit and input-file checksums, and `GEN_*/joblog.txt` has each run's exit status. Set `REEIM_SEED` (default 1) to change the seed stream.
- Runs that fail or are cut off are deleted and replaced, so every generation ends with 2,500 finished runs.
- **If the job hits its time limit, submit exactly the same command again.** It continues from the last finished generation. It first deletes the runs that were cut off, then fills up the generation in progress.
- Do not copy an old `REEIM/out_REEIM/` folder in. The rerun must start from scratch in `out_perch`.

When step 1 is done:

```
sbatch -A <project> REEIM/slurm/pps.sbatch
```

This step reads `REEIM/out_perch/REEIM_res.rds` and writes `REEIM/scenarios_perch/`.
If a PPS job is cut off, submit the same command again: finished scenarios carry a `.done` marker and are skipped, and an unfinished one is rerun.

## 5. Follow progress

```
squeue -u $USER
tail slurm-reeim-abc-<jobid>.out                            # "N out of 2500 particles are done"
for g in REEIM/out_perch/GEN_*; do echo $g $(ls $g | grep -c ^traitStats_); done
```

## 6. Send back

- `REEIM/out_perch/REEIM_res.rds` (and later `REEIM/out_plaice/REEIM_res.rds`).
- The `slurm-*.out` logs.
- For the PPS, the whole `REEIM/scenarios_perch/` folder as a tar file. Check its size first (`du -sh`).

The `GEN_*` folders hold the raw SLiM output. Keep them on the cluster until Petri confirms the results, but they don't need to be sent.

## Files in this folder

| File | What it does |
|---|---|
| `environment.yml` | conda environment (SLiM 5.0, R packages, GNU parallel) |
| `env.sh` | activates it inside jobs; the only file you may need to edit |
| `abc.sbatch` | ABC-SMC for one system (`perch`, `plaice`); resumes after a time limit |
| `ppc.sbatch` | Posterior predictive check of L_T change for a plaice ABC (genetic L_T and survey-estimated PMRN/ogive); resumes |
| `pps.sbatch` | perch posterior predictive simulations |
| `clean_partial.sh` | removes runs cut off by a time limit (called by `abc.sbatch`) |
| `test_one.sh` | single test run |
