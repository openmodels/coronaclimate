## setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

library(dplyr)
library(ggplot2)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

resultdir <- "results"
in.version <- "1201"
out.version <- "1201"
code.version <- "1218"

df <- read.csv("../../cases/panel_all.csv")
df$regid <- paste(df$Country, df$Region, df$Locality)

df2 <- df[, c('regid', 'Country', 'Region', 'Locality', 'population', 'lowest_level', 'implausible')] %>% group_by(regid) %>% summarize(Country=Country[1], Region=Region[1], Locality=Locality[1], population=mean(population), nobs=length(population), lowest_level=min(lowest_level), implausible=max(implausible))
df2.mobile <- df %>% group_by(regid) %>% summarize(portion=mean(!is.na(mobility_pca1[as.Date(Date) > "2020-02-28"])))

for (weight in c('nobs', 'pop', 'region')) {
    for (mobileonly in c(T, F)) {
        results <- read.csv(paste0("../../", resultdir, "/epimodel-", in.version, ".csv"))
	if (mobileonly)
            results <- subset(results, regid %in% df2.mobile$regid[df2.mobile$portion > .5])

        subsuffix <- ifelse(mobileonly, "-mobile", "-all")
        outfile <- paste0("../../", resultdir, "/epimodel-meta-", out.version, subsuffix, "-", weight, ".csv")

        source(paste0("meta-modellib-", code.version, ".R"))
    }
}

