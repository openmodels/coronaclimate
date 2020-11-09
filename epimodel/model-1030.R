## setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

source("../configs.R")

version <- "1030"

casespath <- "../../cases/panel_all.csv"
weather <- c('absh', 't2m', 'tp', 'ssrd', 'utci')
regfilter <- function(rows) T
do.multiproc <- T
mobileonly <- F # Subset it when doing meta analysis

outpath <- paste0("../../results/epimodel-", version, ifelse(mobileonly, "-mobile", ""), ".csv")

source(paste0("modellib-", version, ".R"))


