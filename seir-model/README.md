The files in this directory allow one to fit, monitor, project, and
analyze the SEIR model. Each script is explained in order for
reproducing the results in the paper.

Required R packages:
`dplyr`, `ggplot2`, `rstan`, `lubridate`, `reshape2`, `scales`,
`PBSmapping`, `xtable`, `lfe`, `stringr`, `Hmisc`.

Organization:
Typically, the first line of the script calls `setwd` to this
directory. This path should be updated accordingly.

The following sibling directories to the repository (so, if this
folder is `.../repro/coronavirus/seir-model`, these should be placed
in the `repro` directory) should be created and populated:
 - `cases/`: Put `panel_all.csv` from the Zenodo repository here.
 - `figures`: Empty until figures are generated.
 - `results`: Empty or containing the `results` directory from the Zenodo repository.

Generating the results

 - `model-0314.R` fits the model for each region in the data. It may
   be run individually, or in parallel by calling `multiproc.sh
   model-0314.R [N]`, where `[N]` is a number of processors, on a
   Linux-based system. (`modellib-0314.R` is a helper for this script.)
   
 - Check the progress of this (particularly useful when run in
   parallel) with `check-rhat.R`, `check-params.R`, or
   `check-status.R`.
   
 - Combine the results into a single file with `combine-parallel.R`.
 
 - In cases where a region is attempted, but no initial parameter
   value can be found for the Bayesian MCMC, a placeholder file is
   left. To re-run these regions, these files need to be removed. This
   is performed by `parallel-clean.sh` on a Linux system.
 
 - When parallel processing multiple servers are used, producing
   multiple runs for some regions, use `combine-multidraw.R` after
   producing the single file version.
   
 - `meta-model-0314.R` combines the results from multiple regions
   using the meta-analysis method. (`meta-modellib-0314.R` is a helper
   for this script.)

 - `global-0314.R` fits the single-region model using globally
   aggregated data, as a comparison.

Result analysis and figures:

 - Pairwise comparison (e.g., parameters with and without weather
   estimation): `check-pairwise.R` collects paired results and
   `plot-pairwise.R` displays comparisons.
   
 - `compare-models.R` compares parameter estimates across all results
   from various versions of the model.
   
 - `f-tests.R` calculates significance tests for figure 4.
 
 - Forward projection: `forward-0314.R` provides a function that
   simulates the SEIR model given parameters. `forward-0314.R` is used
   by `projection-impulses.R` to estimate impulse responses; by
   `text-values.R` to produce some values reported in the text; and by
   `validate.R` to check the predictive capacity of the model.  Since
   feedbacks can cause this to diverge from observed results,
   `forward-0314-adaptive.R` provides a corresponding function that
   adjusts case imports to maintain approximate parity. This is used
   by `projection-cities-0314.R` to produce counterfactual case loads
   for cities in `major_cities_selection.csv`.

 - `meta-display.R` produces violin plots of parameter values and the
   mapped effect of a standard deviation change in weather.

 - `plot-climatevar.R` produces figure 5 in the main text, describing
   variance explained by various sources.
   
 - `plot-dynamics.R` plots the region-specific baseline rates, which
   follow a random-walk Bayesian model.

 - `report-noprior.R` creates main figure 3a and `report-model.R` creates main figure 3b.

 - `report-model2.R` creates tables of the model parameters.
