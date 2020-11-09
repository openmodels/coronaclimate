## setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

library(dplyr)
library(ggplot2)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

version <- "1030"

df <- read.csv("../../cases/panel_all.csv")
df$regid <- paste(df$Country, df$Region, df$Locality)

df2 <- df[, c('regid', 'Country', 'Region', 'Locality', 'population', 'lowest_level', 'implausible')] %>% group_by(regid) %>% summarize(Country=Country[1], Region=Region[1], Locality=Locality[1], population=mean(population), nobs=length(population), lowest_level=min(lowest_level), implausible=max(implausible))
df2.mobile <- df %>% group_by(regid) %>% summarize(portion=mean(!is.na(mobility_pca1[as.Date(Date) > "2020-02-28"])))

for (weight in c('pop', 'region', 'nobs')) {
    for (mobileonly in c(F, T)) {
        results <- read.csv(paste0("../../results/epimodel-", version, ".csv"))
	if (mobileonly)
            results <- subset(results, regid %in% df2.mobile$regid[df2.mobile$portion > .5])

        subsuffix <- ifelse(mobileonly, "-mobile", "-all")
        outfile <- paste0("../../results/epimodel-meta-", version, subsuffix, "-", weight, ".csv")

        source(paste0("meta-modellib-", version, ".R"))
    }
}

