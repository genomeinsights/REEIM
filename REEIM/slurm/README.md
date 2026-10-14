# Running the REEIM reruns on a SLURM cluster (UPPMAX)

These notes are for running the REEIM ABC-SMC reruns on a cluster that uses the SLURM scheduler.
Everything runs from one clone of the repository.
Each job uses the cores of one node, and GNU parallel spreads the SLiM runs over them.

## What is run, and in what order

| Step | Job script | Needs | Approx. size |
|---|---|---|---|
| 1. Perch ABC-SMC | `abc.sbatch perch` | nothing | 5 generations × 2,500 SLiM runs ≈ 500 core-hours |
| 2. Perch posterior predictive simulations | `pps.sbatch` | step 1 finished | 24 scenarios × 240 runs ≈ 300–400 core-hours |
| 3. Plaice ABC-SMC | `abc.sbatch plaice` | go-ahead from Petri (inputs still being tested) | similar to step 1 |
| 4. Recovery test | (added later) | step 1; the scripts are rebuilt after it | about 30,000 runs |

The sizes assume about 140 s per SLiM run, which is what an Apple-silicon Mac mini takes.
Cluster cores may be slower, so check with `test_one.sh` first.
Steps 1 and 3 can run at the same time, because each system writes to its own folder.

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

The run takes a few minutes. Expect the `real` time to be roughly 2–3 minutes.
The last line printed should start with 10126, the final year.
Send Petri the time, because it sets the time limit to request.

## 4. Submit

```
sbatch -A <project> -p <partition> -c 64 -t 2-00:00:00 REEIM/slurm/abc.sbatch perch
```

- `-c` is the number of cores. Use a whole node if you can.
- The job writes `REEIM/out_perch/GEN_1 … GEN_5` and, after each generation, `REEIM/out_perch/REEIM_res.rds`.
- **If the job hits its time limit, submit exactly the same command again.** It continues from the last finished generation. It first deletes the runs that were cut off, then fills up the generation in progress.
- Do not copy an old `REEIM/out_REEIM/` folder in. The rerun must start from scratch in `out_perch`.

When step 1 is done:

```
sbatch -A <project> -p <partition> -c 64 REEIM/slurm/pps.sbatch
```

This step reads `REEIM/out_perch/REEIM_res.rds` and writes `REEIM/scenarios_perch/`.
If a PPS job is cut off, delete the newest scenario folder before you resubmit, because scenario folders that already exist are skipped.

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
| `pps.sbatch` | perch posterior predictive simulations |
| `clean_partial.sh` | removes runs cut off by a time limit (called by `abc.sbatch`) |
| `test_one.sh` | single test run |
