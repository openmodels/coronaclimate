cp ../../results/* ../../results-saved
Rscript combine-parallel.R
Rscript meta-model-0314.R

## Make figure 2
Rscript projection-impulses.R
Rscript report-model.R
