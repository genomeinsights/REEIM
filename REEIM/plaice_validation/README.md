# North Sea plaice hold-out validation

A second study system for REEIM, used as a validation. The model gets only contemporary data
(2015–2020) plus the relative exploitation history. Its output is then checked against two
quantities it never sees:

- **(a)** the assessed harvest rate H, catch / total biomass in 2020 = 0.078 (ICES 2025 assessment), and
- **(b)** the published decline in female maturation size (PMRN Lp50), e.g. age 4: 34.4 → 30.1 cm
  (cohorts 1955–64 → 1985–94, Grift et al. 2003) and 41.6 → 22.9 cm (1900s → 2000s, van Walraven et al. 2010).

Background, choices and scientific caveats are in the project note `plaice_holdout_feasibility.md`.
Issue numbers (I-xx) refer to the project's `issue_register.md`.

## Layout

| Path | Contents |
|---|---|
| `../input/plaice/` | Model inputs in the standard format (see `../input/README.md`), plus `growth_priors.csv` for the growth fit. **Main set:** g fixed at 0.25. |
| `../input/plaice_Hfree/` | Identical to `../input/plaice/` except the H prior: 0–0.2 instead of 0.06–0.11 (ICES catch/TSB). Run alongside the main set to show what the data say about H without the assessment. |
| `../input/plaice_gfree/` | Growth constants from the free-g fit (g ≈ 0.81) (I-18). Kept for reference only and not run: g is set at a realistic 0.25 because size at age does not identify it. |
| `data_prep/get_observed_data_plaice.R` | DATRAS BTS Q3 and NS-IBTS Q1, 2015–2020 → `data/observed_data_plaice_{female,pooled}.RData`, maturity ogive, summary CSVs. |
| `data_prep/build_fishing_intensity.py` | Relative fishing intensity 1850–2020 → `../input/plaice/fishing_intensity.csv`. |
| `data_prep/sources/` | ICES 2024 assessment summary, Rijnsdorp & Millner (1996) Table 1, and the derived pre-1957 F on the ICES scale. |
| `growth_fit/deterministic_ABC.R` | Copy of `deterministic_ABC/deterministic_ABC.R` that reads L_T, max age and priors from a system input folder (I-19). With no input folder it behaves exactly like the perch script. |
| `growth_fit/output/` | Plaice growth-fit results (2,500 particles × 5 generations). |
| `validation_targets/` | Published Lp50 values and the ICES 2025 assessment series with catch / TSB. |
| `data/` | Generated observed data, and the cached raw DATRAS downloads (`raw_*.rds`). |
| `ppc_plaice.R` | Posterior predictive check of evolutionary change after a plaice ABC (run by `../slurm/ppc.sbatch`; see below). |
| `test_runs/` | Pre-ABC test runs and their scripts (outputs not in git): prior draws (`make_test_runs.py`), S_LT × w grid (`make_grid.py`, `grid_LT_change.R`), maximum-response stress test (`make_stress.py`, `analyse_stress.R`), virtual-ecologist survey model (`REEIM_vesurvey.slim`, `make_vesurvey.py`, `analyse_vesurvey.R`), plots against fishing intensity (`plot_LT_vs_fishing.R`, `plot_timing_signature.R`), and the ABC-pilot diagnostics (`abc_pilot/size_diag.R`, `size_fix.R`, `compat.R`). |

## Running it (from the repository root)

```sh
# 1. Observed data. Uses the cached raw_*.rds files in data/ if present; otherwise needs
#    internet access and install.packages("icesDatras").
Rscript REEIM/plaice_validation/data_prep/get_observed_data_plaice.R

# 2. Fishing trajectory (optional; already in ../input/plaice). Arguments: war-year fraction, terminal year.
python3 REEIM/plaice_validation/data_prep/build_fishing_intensity.py \
  REEIM/plaice_validation/data_prep/sources REEIM/input/plaice/fishing_intensity.csv 0.1 2020

# 3. Deterministic growth fit (optional; results already in growth_fit/output)
Rscript REEIM/plaice_validation/growth_fit/deterministic_ABC.R \
  REEIM/plaice_validation/data/observed_data_plaice_female.RData 4 \
  REEIM/input/plaice REEIM/plaice_validation/growth_fit/output

# 4. REEIM ABC-SMC (copy REEIM/input/ along with the model when running on another machine)
Rscript REEIM/run_REEIM_ABC.R . <cores> ./REEIM/input/plaice
Rscript REEIM/run_REEIM_ABC.R . <cores> ./REEIM/input/plaice_Hfree     # same, H prior 0-0.2
```

## Data (DATRAS, 2015–2020)

- **Size structure and size at age:** BTS quarter 3 beam-trawl survey. North Sea only: hauls at 51° N or below (the Channel) and west of Great Britain (the Irish Sea and Bristol Channel, a different stock) are dropped (I-20). That leaves 1,702 valid hauls and 148,505 plaice measured. Catches are standardised per metre of beam per hour, because the vessels use 4, 7 and 8 m beams. The otoliths are length-stratified, so size at age comes from an age–length key applied to catch at length.
- **Females:** the length frequencies are mostly unsexed, so the female size distribution is catch at length × the female share at that length in the otolith samples. Female size at age is 107, 163, 211, 247, 274 and 301 mm at ages 0–5, levelling off near 350–360 mm.
- **Maturity:** NS-IBTS quarter 1 (plaice spawn in Q1), North Sea and Skagerrak with the Kattegat excluded. Female L50 = 219 mm (L25–L75 170–269, n = 5,792); male L50 = 156 mm. `LT_obs` = 219.
- **Size-distribution window:** 12–50 cm (`kde_from` = 12), because the beam trawl under-catches smaller age-0 plaice (I-07).

## Growth constants (I-18, I-19)

The deterministic growth ABC was run on the female size at age, with L_T = 219 mm.

| Model | h (mm/yr) | h0 | g_i | g | final ε |
|---|---|---|---|---|---|
| (5) L_T fixed, g free | 99.4 (58.6–143.4) | 1.16 | 0.73 | 0.81 (0.36–1.32) | 0.503 |
| (6) L_T fixed, g = 0.25 | 50.0 (42.7–57.1) | 2.36 | 0.11 | 0.25 | 0.484 |

g is only weakly identified from size at age: fixing it at a biologically plausible 0.25 fits as
well as letting it go free (least-squares SRMSE 0.437 vs 0.426). The free fit pushes the growth
slow-down at the record 2010s stock size into g. The two curves differ mainly beyond age 8, where
the data are sparse. The main set therefore uses g = 0.25, and the free-g set is kept for
reference only (not run). `hmaxPrior` = h × 1.616 (1 + γ·D_K^δ), as for perch (38.3 → 62).

Priors that depend on trait scale (S_hmax, S_LT, f) are the perch ranges rescaled by the ratio of
the plaice to the perch reference values. They should be checked on test runs.

## Size-distribution window (kde_from = 16 cm)

BTS Q3 is an offshore beam-trawl survey and under-catches plaice below about 16 cm (age 0–1 still in
coastal nurseries, I-07). With the comparison starting at 12 cm, 58% of the squared size-distribution
error of the best pilot particles came from 12–16 cm, while size at age 1–5 and the age composition
fitted well. Starting at 16 cm brings the size-distribution SRMSE of those particles from 0.39 to 0.12,
the level of perch (`test_runs/abc_pilot/size_diag.R`, `size_fix.R`).

## Evolutionary change: posterior predictive check

The published declines are PMRN Lp50 estimates from survey and market samples, not the genetic L_T, so
the check is made on both:

1. **Genetic L_T** (`ppc_plaice.R genetic`): weighted change over the accepted particles of the last ABC
   generation, from their `traitStats` files; no new runs.
2. **Survey statistics** (`ppc_plaice.R prepare` / `survey`): 250 posterior draws are rerun with
   `test_runs/REEIM_vesurvey.slim`, a copy of `REEIM.slim` that also writes a random sample of 2,000
   females (age ≥ 1) every year 1900–2020 (age, length, mature, first-time maturing; about 4 MB per run).
   To each simulated survey are fitted the maturity ogive L50 (per year) and the PMRN Lp50 at age 4
   (per decade), the latter both with first-time maturing known and with the Barot et al. (2004)
   demographic method from age-specific ogives and growth, as in Grift et al. (2003) and van Walraven
   et al. (2010).

3. **Evolutionary sensitivity** (`ppc_sensitivity_summary.csv`): each posterior draw is also run with
   H = 0 and the same seed (identical burn-in). The difference in the genetic L_T change between the pair
   is the effect of fishing; divided by 10 H it gives the fishing-induced change in L_T in % per 0.1 of
   terminal harvest rate. This separates how hard the stock was fished (H, from the assessment) from how
   strongly this population responds (S_LT, f, w, Shmax, constrained by the present-day data), which trade
   off against H when H is not known. `ppc_plaice.R genetic` also reports a naive version (change / 10 H,
   assuming no change without fishing) from the accepted particles alone.

For each statistic, `ppc_*_summary.csv` gives the median, 90% interval and the posterior predictive
probability of a decline at least as large as published, for 1905–85 (−20%), 1960–90 (−12.5%),
1985–2005 (−31%) and 1900–2020; `fig/ppc_trajectories.png` shows the bands with fishing intensity.
Keep `REEIM_vesurvey.slim` in step with `REEIM.slim`: it differs only in the three blocks marked
"VE survey".

Supporting evidence from the test runs (`test_runs/abc_pilot/compat.R`): among 598 runs scored with the
16 cm window, the best fit reachable falls sharply with the size of the simulated 1960–90 decline
(total distance 0.59 with no decline, 1.26 at ≥ 4%, 3.23 at ≥ 12.5%), and all three distance
components get worse together.
