setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

outfile <- "../../results/epimodel-meta-1217-all-nobs.csv"
outfile.mob <- "../../results/epimodel-meta-1217-mobile-nobs.csv"

allrecorded <- read.csv(outfile)
results <- subset(allrecorded, Country != "" & group == "Combined" & Region == "" & Locality == "")

library(dplyr)

results2 <- results %>% group_by(regid) %>% summarize(delay=mu[param == 'invkappa'] + mu[param == 'invsigma'])
