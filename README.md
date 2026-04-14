# Retrospective Eco–Evolutionary Inference Model (REEIM)

This repository contains all code and data required to reproduce the analyses in *Inferring human-induced evolution where history is missing* which uses a retrospective eco–evolutionary inference framework (REEIM) to infer historical fishing intensity and associated life-history evolution from contemporary monitoring data.

The workflow integrates empirical data processing, eco-genetic individual-based simulations in **SLiM**, Approximate Bayesian Computation with Sequential Monte Carlo (ABC–SMC), and posterior predictive simulations (PPS).

------------------------------------------------------------------------

## Repository structure

### `./empirical_data/`

Contains all empirical inputs used for model calibration.

-   **`R/get_observed_data.R`**\
    Processes raw monitoring data and extracts the summary statistics used for ABC inference.\
    Output:

    -   `./empirical_data/observed_data.RData`\
        → the *contemporary life-history patterns* (size structure, growth, and maturation) used as ABC targets.

These data represent the present-day demographic patterns the model is calibrated against.

------------------------------------------------------------------------

### `./R/SSF.R`

Prepares fishing selectivity and historical fishing intensity inputs for the SLiM simulations, using empirical catch and vulnerability data stored in `./empirical_data/`.

These inputs define the **fixed fishing regime**:

-   size-dependent vulnerability to capture

-   historical development of fishing intensity

These components are treated as fixed (non-inferred) in all simulations.

------------------------------------------------------------------------

### `./REEIM/`

Core simulation and inference pipeline.

-   **`REEIM.slim`**\
    Eco-genetic individual-based model implemented in SLiM.\
    Used for:

    -   ABC–SMC inference

    -   posterior predictive simulations

-   **`run_REEIM_ABC.R`**\
    Runs the full ABC–SMC pipeline:

    -   launches SLiM simulations

    -   computes distances to observed data

    -   stores accepted particles per generation

    Output:

    -   `./REEIM/out_REEIM/` — simulation output for each ABC generation

    -   `./REEIM/REEIM_res.rds` — consolidated object containing all particles, weights, and diagnostics

-   **`RUN_PPS.R`**\
    Uses `REEIM_res.rds` to run posterior predictive simulations (PPS), conditioning on selected harvest intensities and propagating demographic and genetic stochasticity.

------------------------------------------------------------------------

### `./R/`

Analysis and figure-generation scripts.

-   **`diagnostics.R`**\
    ABC–SMC diagnostics (tolerance reduction, ESS, acceptance rates).

-   **`PPC_ABC_PPS.R`**\
    Posterior predictive checks (PPCs) for the first five ABC generations as well as for the posterior predictive simulations (PPS)

-   **`PPS.R`**\
    Analysis of posterior predictive simulations, including scenario comparisons.

-   **`PPS_haldanes.R`**\
    Analysis of posterior predictive simulations that includes transformation of evolutionary change to Haldanes.

-   **`posterior_contraction.R`**\
    Visualization of prior–posterior contraction across ABC generations.

-   **`Table_ABC_params.R`**\
    Generates posterior summaries of inferred parameters (Table S2).

-   **`trait_cor_with_H.R`**\
    Quantifies and visualizes evolutionary trait responses as a function of inferred harvest intensity (H).

------------------------------------------------------------------------

## Conceptual overview

The REEIM framework proceeds in three stages:

1.  **Empirical calibration**\
    Contemporary life-history patterns are extracted from monitoring data and used as ABC targets.

2.  **Retrospective inference (ABC–SMC)**\
    Historical harvest intensity and eco-genetic dynamics are inferred by fitting SLiM simulations to observed patterns.

3.  **Posterior predictive simulations (PPS)**\
    Forward simulations conditional on the inferred posterior are used to validate model realism and explore counterfactual scenarios.

------------------------------------------------------------------------

## Reproducibility notes

-   All stochastic simulations are fully parameterized via saved particle files.

-   ABC–SMC results are fully reproducible from `REEIM_res.rds`.

-   Posterior predictive simulations reuse the same SLiM model and differ only in parameter conditioning.

------------------------------------------------------------------------

## Contact

For questions regarding the code or analyses, please contact:\
[**petri\@genomeinsights.fi**](mailto:petri@genomeinsights.fi)
