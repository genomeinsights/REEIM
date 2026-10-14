# Study-system input files

Everything that describes a particular study system lives in one folder here
(`perch/` is the system analysed in the paper). `REEIM.slim` reads the folder
given with `-d inputDir="..."`, and `run_REEIM_ABC.R` and `RUN_PPS.R` pass it
for you. Without `-d inputDir`, SLiM uses `./REEIM/input/perch` relative to
the folder it is started from, so run it from the repository root (as all the
R drivers do) or copy `REEIM/input/` along with `REEIM/REEIM.slim`. To set up a new system, copy `perch/` to a new folder, edit the five
files and run

```
Rscript REEIM/run_REEIM_ABC.R <repo folder> <cores> ./REEIM/input/<system>
```

| File | Read by | Contents |
|---|---|---|
| `parameters.csv` | `REEIM.slim` | Species and system constants: carrying capacity, maximum age, reference trait values, growth and density dependence, genetics, length–weight relationship, burn-in length, burn-in abort limits and how the selectivity curve is extended. Any value given on the command line with `-d` overrides the file. |
| `fishing_intensity.csv` | `REEIM.slim`, `run_REEIM_ABC.R` | Relative fishing intensity per year (`year,intensity`), scaled so the maximum is 1. Fishing starts after `Tdelta` burn-in cycles. The last row is the present-day observation year, and `H` is the harvest rate in that year. The first row is used only to fix the timing: the year after burn-in uses the second row. |
| `selectivity.csv` | `REEIM.slim` | Size-selectivity curve (`size_class,vulnerability`). A fish of length L mm is in class `trunc(L / selBinWidth)`, so with 10 mm bins class 6 is 60–69 mm. Classes must be consecutive. Smaller fish get `selBelow` and larger fish get `selAbove`; set `selAbove` to 0 for a dome-shaped (e.g. gillnet) curve. |
| `priors.csv` | `run_REEIM_ABC.R` | Uniform priors for the inferred parameters (`parameter,lower,upper`). |
| `abc_settings.csv` | `run_REEIM_ABC.R` | Observed data file, observed mean size at maturation `LT_obs`, and the bandwidth and range (cm) of the size-distribution comparison. |

`parameters.csv` and `abc_settings.csv` are `name,value,description` tables;
lines starting with `#` are comments. In `parameters.csv`, write integers
without a decimal point (`10000`) and reals with one (`62.0`), as SLiM does
for `-d` values.

The ABC compares simulations with the observed data at cycle
`Tdelta + nrow(fishing_intensity.csv) - 1` (10125 for perch), which
`run_REEIM_ABC.R` computes from these files.

## Perch values

`fishing_intensity.csv` and `selectivity.csv` are written by `R/SSF.R`
from the angler survey and the vulnerability data in `empirical_data/`. The
values in `parameters.csv` are the ones used for the ABC–SMC in the paper.
`gi` and `h0` come from the deterministic Lester model fit, and the
length–weight coefficients are those passed by `run_REEIM_ABC.R`
(a = 3.95801e-6, b = 3.18606). The posterior predictive simulations used the
earlier `REEIM.slim` defaults (a = 3.958007e-6, b = 3.186062); `RUN_PPS.R`
still passes those so the published simulations are reproduced.
