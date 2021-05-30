setwd("~/research/coronavirus/code/epimodel")

library(dplyr)
library(reshape2)

weathervars <- c('t2m', 'tp', 'ssrd', 'utci')
version <- "0314"

## Let's do everything by bootstrapping!
## Compare it always to some plausible counterfactual.

if (F) {
    ## Generate fake data for testing
    version <- 'test' #"0314"

    make.row <- function(regid, param, dist) {
        data.frame(regid, param, mu=mean(dist), ci2.5=quantile(dist, .025), ci25=quantile(dist, .25),
                   ci50=quantile(dist, .5), ci75=quantile(dist, .75), ci97.5=quantile(dist, .975),
                   group='Raw')
    }

    do.null <- F

    for (do.weather in c(F, T)) {
        results <- data.frame()
        for (regid in paste0("Nevernever ", 1:100)) {
            if (do.null)
                results <- rbind(results, make.row(regid, 'error', abs(rnorm(1000, runif(1), runif(1)))))
            else
                results <- rbind(results, make.row(regid, 'error', abs(rnorm(1000, runif(1 + !do.weather), runif(1)))))
            if (do.weather)
                for (param in c(paste0('e.', weathervars), paste0('o.', weathervars)))
                    if (do.null)
                        results <- rbind(results, make.row(regid, param, rnorm(1000, 0, .1 * runif(1))))
                    else
                        results <- rbind(results, make.row(regid, param, rnorm(1000, rnorm(1, 0, .1), .1 * runif(1))))
            if (do.null)
                dynamics <- data.frame(logbeta=cumsum(rnorm(365)), logomega=c(NA, cumsum(rnorm(364))))
            else
                dynamics <- data.frame(logbeta=cumsum(rnorm(365)) + sin(2 * pi * (1:365) / 365), logomega=c(NA, cumsum(rnorm(364))) + sin(2 * pi * (1:365) / 365))
            write.csv(dynamics, paste0("../../results/epimodel-meta-", version, "-noprior-dynamics.csv-", regid), row.names=F)
        }

        if (do.weather)
            write.csv(results, paste0("../../results/epimodel-meta-", version, "-noprior-all-nobs-nodel.csv"), row.names=F)
        else
            write.csv(results, paste0("../../results/epimodel-meta-", version, "-noweather-all-nobs-nodel.csv"), row.names=F)
    }
}

## 1. Existence of any weather effect
## Compare noweather to noprior: Is error less?

## Version 1
results.np <- subset(read.csv(paste0("../../results/epimodel-meta-", version, "-noprior-all-nobs-nodel.csv")), group == "Raw" & param == 'error')
results.nw <- subset(read.csv(paste0("../../results/epimodel-meta-", version, "-noweather-all-nobs-nodel.csv")), group == "Raw" & param == 'error')

draws.np <- results.np[, c('regid', "ci2.5","ci25","ci50","ci75","ci97.5")]
draws.nw <- results.nw[, c('regid', "ci2.5","ci25","ci50","ci75","ci97.5")]

draws <- draws.np %>% left_join(draws.nw, by='regid', suffix=c('.np', '.nw'))
draws.np2 <- melt(draws[, c('regid', paste0(c("ci2.5","ci25","ci50","ci75","ci97.5"), '.np'))], id='regid')
draws.nw2 <- melt(draws[, c('regid', paste0(c("ci2.5","ci25","ci50","ci75","ci97.5"), '.nw'))], id='regid')

draws.np2$compare <- draws.nw2$value

pval.weather1 <- mean(draws.np2$value > draws.np2$compare, na.rm=T)

draws.np2$rss2 <- draws.np2$value^2 * 365
draws.np2$rss1 <- draws.np2$compare^2 * 365

df1 <- 8
df2 <- 2*365 - (10 + 8 + 2*6 + 1 + 365)
fstat <- ((draws.np2$rss1 - draws.np2$rss2) / df1) / (draws.np2$rss2 / df2)
pval.weather2 <- mean(1 - pf(fstat, df1, df2), na.rm=T)

## Version 2
results.np <- subset(read.csv(paste0("../../results/epimodel-meta-", version, "-noprior-all-nobs-nodel.csv")), group == "Combined" & regid == '  ' & param == 'error')
results.nw <- subset(read.csv(paste0("../../results/epimodel-meta-", version, "-noweather-all-nobs-nodel.csv")), group == "Combined" & regid == '  ' &  param == 'error')

## Error considerably less with no weather

## 2. Existence of a transmission/detection channel
## Null hypothesis: same distribution, but each mean 0.

do.set <- "Combined"

for (prefix in c('e.', 'o.')) {
    results <- subset(read.csv(paste0("../../results/epimodel-meta-", version, "-noprior-all-nobs-nodel.csv")), group == do.set & param %in% paste0(prefix, weathervars))

    alldraws <- c()
    for (regid in unique(results$regid)) {
        subres <- results[results$regid == regid,]
        draws <- subres[, c("ci2.5","ci25","ci50","ci75","ci97.5"),]

        ## Under null hypothesis p(all < 0) is .5^K, so multiply by 2^K

        ## Assume no correlation
        alldraws <- c(alldraws, sapply(1:1000, function(ii) all(sign(sapply(1:nrow(draws), function(ii) draws[ii, sample.int(5, 1)])) != sign(subres$mu))))
    }

    if (prefix == 'e.')
        pval.weather.trans <- (2^nrow(draws)) * mean(alldraws) # test passed
    else
        pval.weather.detect <- (2^nrow(draws)) * mean(alldraws) # test passed
}

pval.weather.trans
pval.weather.detect

## 3. Behaviour affects transmission/detection
## Null hypothesis: logomega and logbeta are random walk.
results <- read.csv(paste0("../../results/epimodel-", version, "-noprior-nodel.csv"))
regids <- unique(results$regid)

logbeta.pvals <- c()
logomega.pvals <- c()
for (regid in regids) {
    dynpath <- paste0("../../results/epimodel-", version, "-noprior-dynamics.csv-", regid)
    if (!file.exists(dynpath))
        next

    dynamics <- read.csv(dynpath)

    for (col in c('logbeta', 'logomega')) {
        dyndiff <- diff(dynamics[, col])
        pval <- cor.test(dyndiff[-1], dyndiff[-length(dyndiff)], alternative='t', na.rm=T)$p.value

        if (col == 'logbeta')
            logbeta.pvals <- c(logbeta.pvals, pval)
        else
            logomega.pvals <- c(logomega.pvals, pval)
    }
}

pval.behave.trans <- mean(logbeta.pvals, na.rm=T) # failed in non-null
pval.behave.detect <- mean(logomega.pvals, na.rm=T) # failed in non-null

library(ggplot2)

ggplot(data.frame(pval=c(pval.behave.trans, pval.behave.detect),
                  var=rep(c('log(beta)', 'log(gamma)'), each=length(pval.behave.trans))),
       aes(pval)) +
    facet_wrap(~ var) + geom_histogram() +
    scale_x_log10()

## 4. Governance/climate affects parameters
## Null hypothesis: Parameters not explained by governance/climate

## Calculated in plot-climatevar.R
