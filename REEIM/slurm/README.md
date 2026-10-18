# Running the REEIM reruns on a SLURM cluster (Dardel, PDC)

These notes are for running the REEIM ABC-SMC reruns on a cluster that uses the SLURM scheduler.
Everything runs from one clone of the repository.
Each job uses the cores of one node, and GNU parallel spreads the SLiM runs over them.

## What is run, and in what order

| Step | Job script | Needs | Approx. size |
|---|---|---|---|
| 1. Perch ABC-SMC | `abc.sbatch perch` | nothing | 5 generations × 2,500 SLiM runs (≈ 14,000 runs incl. replaced ones) ≈ 2,800 core-hours, about 22 node-hours: plan on two 24 h jobs |
| 2. Perch posterior predictive simulations | `pps.sbatch` | step 1 finished | 24 scenarios × 240 runs ≈ 1,000 core-hours (if the runs are as long as the ABC runs) |
| 3. Plaice ABC-SMC (two jobs) | `abc.sbatch plaice` (H prior 0.06–0.11, ICES catch/TSB) and `abc.sbatch plaice_Hfree` (H prior 0–0.2) | go-ahead from Petri (inputs still being tested) | about 2,800 core-hours each if the runs are as long as perch's |
| 3b. Plaice posterior predictive check | `ppc.sbatch plaice` and `ppc.sbatch plaice_Hfree` | the matching step-3 job finished | 250 draws + 250 paired H = 0 controls = 500 runs each ≈ 90 core-hours, about 1 node-hour |
| 3c. Plaice evolutionary-change scenarios | `scenarios.sbatch plaice` | the matching step-3 job finished | 50 draws × 5 genetics × 5 H = 1,250 runs ≈ 250 core-hours, about 2–3 node-hours (high-mutation runs are slower) |
| 4. Recovery test | (added later) | step 1; the scripts are rebuilt after it | about 30,000 runs ≈ 5,000 core-hours |

The sizes are measured on Dardel (October 2026, one `main` node, the `test_one.sh` parameter set): with 256 runs at once (two per physical core) a run takes about 17 min and a node finishes about 760 runs per hour, i.e. about 0.17 core-hours per run. With 128 runs at once (one per core) a run takes about 12 min but a node finishes only about 550 runs per hour. The estimates add about 15% for idle cores at the end of each batch of runs (since the batch change, one batch per generation). Run times vary with the parameter values, so treat them as ±30%.
Check with `test_one.sh` first.

**Dardel specifics.** A thin node has 128 physical cores (256 logical CPUs) and 256 GB of memory.
The `main` partition gives whole nodes, is charged for all 128 cores, and allows at most 24 hours per job.
The job scripts start as many SLiM runs at once as `-c` (`SLURM_CPUS_PER_TASK`). On Dardel `-c` counts hardware threads, two per physical core, and a `main` job gets the whole node either way. The job scripts therefore default to `-p main -c 256 -t 24:00:00`: 256 SLiM runs at once, which gives about 40% more runs per node-hour than 128 at once (`-c 128`).
With that, one ABC generation (2,500 runs) takes about 3.5–4 hours on one node (GEN_1 longer, because runs with invalid output are replaced), so the five generations need about 22 hours: submit the job together with a `--dependency=afterany` continuation (see below).
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

The run takes a few minutes. Expect the `real` time to be roughly 10–15 minutes on Dardel (3–6 minutes on a recent desktop).
The last line printed should start with 10126, the final year.
Send Petri the time, because it sets the time limit to request.

## 4. Submit

```
sbatch -A <project> REEIM/slurm/abc.sbatch perch
```

- The defaults are one whole `main` node running 256 SLiM runs at once, for 24 hours (see Dardel specifics above). Override the defaults on the command line if needed, e.g. `-p shared -c 32` for a small test.
- To queue a continuation automatically, submit a second job that starts when the first ends: `sbatch -A <project> --dependency=afterany:<jobid> REEIM/slurm/abc.sbatch perch`.
- The job writes `REEIM/out_perch/GEN_1 … GEN_5` and, after each generation, `REEIM/out_perch/REEIM_res.rds`.
- Every SLiM run gets an explicit seed. `REEIM/out_perch/manifest.tsv` lists each run's generation, seed and full command, `manifest_header_*.txt` records the commit and input-file checksums, `GEN_*/batch_<k>.txt` holds each batch's commands and `GEN_*/joblog_<k>.txt` each run's exit status (runs from before the batch change are in `GEN_*/joblog.txt`). A batch cut off by a time limit is finished with the same commands when the job is resubmitted. Set `REEIM_SEED` (default 1) to change the seed stream. For a short test run, set `REEIM_N_PARTICLES` and `REEIM_GENS` (defaults 2500 and 5), e.g. `sbatch -A <project> -p shared -c 32 -t 4:00:00 --export=ALL,REEIM_N_PARTICLES=100,REEIM_GENS=2,REEIM_OUT=./REEIM/test_perch/ REEIM/slurm/abc.sbatch perch`. `REEIM_OUT` keeps the test out of `out_perch`; delete the test folder before the real run.
- Runs that fail or are cut off are deleted and replaced, so every generation ends with 2,500 finished runs.
- **If the job hits its time limit, submit exactly the same command again.** It continues from the last finished generation. It first deletes the runs that were cut off, then fills up the generation in progress.
- Do not copy an old `REEIM/out_REEIM/` folder in. The rerun must start from scratch in `out_perch`.

When step 1 is done:

```
sbatch -A <project> REEIM/slurm/pps.sbatch
```

This step reads `REEIM/out_perch/REEIM_res.rds` and writes `REEIM/scenarios_perch/`.
If a PPS job is cut off, submit the same command again: finished scenarios carry a `.done` marker and are skipped, and an unfinished one is rerun.

Plaice runs the same way, with two ABCs that can run at the same time as each other and as perch:

```
sbatch -A <project> REEIM/slurm/abc.sbatch plaice
sbatch -A <project> REEIM/slurm/abc.sbatch plaice_Hfree
```

When each has finished, run its posterior predictive check, and the scenarios for the main set:

```
sbatch -A <project> REEIM/slurm/ppc.sbatch plaice
sbatch -A <project> REEIM/slurm/ppc.sbatch plaice_Hfree
sbatch -A <project> REEIM/slurm/scenarios.sbatch plaice
```

They read `REEIM/out_<system>/` and write `REEIM/ppc_<system>/` and `REEIM/scen_<system>/`. As for the ABC, every run has an explicit seed (in `cmds.txt`), `manifest_header_*.txt` records the commit and input checksums, and `joblog.txt` each run's exit status; `REEIM_SEED` changes the seed stream, and `REEIM_ABC_OUT` / `REEIM_PPC_OUT` / `REEIM_SCEN_OUT` point a test run at other folders. If a job is cut off, submit the same command again: finished runs are skipped.

## 5. Follow progress

```
squeue -u $USER
tail slurm-reeim-abc-<jobid>.out                            # "N out of 2500 particles are done"
for g in REEIM/out_perch/GEN_*; do echo $g $(ls $g | grep -c ^traitStats_); done
```

## 6. Send back

- `REEIM/out_perch/REEIM_res.rds`, `REEIM/out_plaice/REEIM_res.rds` and `REEIM/out_plaice_Hfree/REEIM_res.rds`.
- The `slurm-*.out` logs.
- For the PPS, the whole `REEIM/scenarios_perch/` folder as a tar file. Check its size first (`du -sh`).
- For plaice, `REEIM/ppc_<system>/` and `REEIM/scen_<system>/` without their `out/` folders (the summaries, `fig/`, `design.csv`, `cmds.txt`, `joblog.txt` and `manifest_header_*.txt`); the `out/` folders hold about 4 MB per run of survey output.

The `GEN_*` folders hold the raw SLiM output. Keep them on the cluster until Petri confirms the results, but they don't need to be sent.

## Files in this folder

| File | What it does |
|---|---|
| `environment.yml` | conda environment (SLiM 5.0, R packages, GNU parallel) |
| `env.sh` | activates it inside jobs; the only file you may need to edit |
| `abc.sbatch` | ABC-SMC for one system (`perch`, `plaice`); resumes after a time limit |
| `ppc.sbatch` | Posterior predictive check of L_T change for a plaice ABC (genetic L_T and survey-estimated PMRN/ogive); resumes |
| `scenarios.sbatch` | Expected change in L_T for posterior draws × additive genetic variance × H, with fit to the present-day data; resumes |
| `pps.sbatch` | perch posterior predictive simulations |
| `clean_partial.sh` | removes runs cut off by a time limit (called by `abc.sbatch`) |
| `test_one.sh` | single test run |
