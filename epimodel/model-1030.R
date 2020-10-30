## setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

source("../configs.R")

version <- "1030"

casespath <- "../../cases/panel_all.csv"
weather <- c('absh', 't2m', 'tp', 'ssrd', 'utci')
regfilter <- function(rows) T
do.multiproc <- T

for (mobileonly in c(F, T)) {
    outpath <- paste0("../../results/epimodel-", version, ifelse(mobileonly, "-mobile", "-all"), ".csv")

    source(paste0("modellib-", version))
}

