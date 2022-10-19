setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

rawres <- read.csv('../../results/epimodel-0907.csv')
results <- read.csv('../../results/epimodel-meta-0907.csv')
subres <- subset(results, regid == "  ")

casespath <- "~/Downloads/panel-prepped_MLI.RData"
load(casespath)

## Estimate SD variations
library(lfe)

weather <- c('absh', 'r', 't2m', 'tp')
dm.weather <- demeanlist(df[, weather], list(factor(df$regid)))

sds <- as.numeric(apply(dm.weather, 2, sd))

logbeta <- mean(rawres$mu[rawres$param == 'logbeta'])
omega <- subres$mu[subres$param == 'omega']

e.mu <- c(subres$mu[subres$param == 'e.absh'], subres$mu[subres$param == 'e.r'],
          subres$mu[subres$param == 'e.t2m'], subres$mu[subres$param == 'e.tp'])
o.mu <- c(subres$mu[subres$param == 'o.absh'], subres$mu[subres$param == 'o.r'],
          subres$mu[subres$param == 'o.t2m'], subres$mu[subres$param == 'o.tp'])

exp(logbeta + e.mu * sds)
omega + o.mu * sds

tbl <- data.frame(weather=c('Baseline', weather), prob.trans=100 * exp(logbeta + c(0, e.mu * sds)),
                  perc.trans=100 * (c(NA, exp(logbeta + e.mu * sds) / exp(logbeta)) - 1),
                  prob.detect=100 * (omega + c(0, o.mu * sds)), perc.detect=100 * ((omega + c(NA, o.mu * sds)) / omega - 1))

library(xtable)
print(xtable(tbl, digits=2), include.rownames=F)
