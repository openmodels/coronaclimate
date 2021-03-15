setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

priordf <- rbind(read.csv("code/epimodel/coeffs_pol1_distlag_c02_1-12_global_modelselection-A00_nomob_10.csv"),
                 read.csv("code/epimodel/coeffs_pol1_distlag_c02_1-12_global_modelselection-A00_nomob_11.csv"))

source("forward-0210.R")

## Average these up over the periods

library(dplyr)
res.avg <- results %>% group_by(delaysum, pulsetime, effect) %>%
    summarize(estimate=c(mean(impresp[time - pulsetime > 1.5 & time - pulsetime < 3.5]),
                          mean(impresp[time - pulsetime > 3.5 & time - pulsetime < 12.5]),
                          mean(impresp[time - pulsetime > 12.5 & time - pulsetime < 30])),
              period=c('one', 'two', 'three'))
res.avg$period <- factor(res.avg$period, levels=c('one', 'two', 'three'))

ggplot(res.avg, aes(delaysum, estimate, colour=factor(pulsetime), shape=effect)) +
    geom_point()

ggplot(subset(res.avg, pulsetime == 100), aes(delaysum, estimate, colour=effect)) +
    facet_wrap(~ period) +
    geom_point()

res.avg$group <- paste(res.avg$effect, res.avg$period, res.avg$pulsetime)
res.avg$logest <- log(res.avg$estimate)
res.avg$logest[!is.finite(res.avg$logest)] <- NA

res.avg$predicted <- NA
transcoeff <- data.frame()
for (group in unique(res.avg$group)) {
    if (all(is.na(res.avg$logest[res.avg$group == group])))
        next
    if (length(grep('three', group)) > 0)
        next
    if (length(grep('250', group)) > 0)
        next
    if (median(sign(res.avg$estimate[res.avg$group == group])) == -1) {
        res.avg$logest[res.avg$group == group] <- log(-res.avg$estimate[res.avg$group == group])
        res.avg$logest[res.avg$group == group & !is.finite(res.avg$logest)] <- NA
    }
    mod <- lm(logest ~ delaysum, data=res.avg[res.avg$group == group,])
    if (median(sign(res.avg$estimate[res.avg$group == group])) == -1) {
        res.avg$predicted[res.avg$group == group] <- -exp(predict(mod, res.avg[res.avg$group == group,]) + var(mod$resid) / 2)
        transcoeff <- rbind(transcoeff, data.frame(group, sign=-1, intercept=mod$coeff[1] + var(mod$resid) / 2, slope=mod$coeff[2]))
    } else {
        res.avg$predicted[res.avg$group == group] <- exp(predict(mod, res.avg[res.avg$group == group,]) + var(mod$resid) / 2)
        transcoeff <- rbind(transcoeff, data.frame(group, sign=1, intercept=mod$coeff[1] + var(mod$resid) / 2, slope=mod$coeff[2]))
    }
}

ggplot(subset(res.avg, pulsetime == 100), aes(delaysum, estimate, colour=effect)) +
    facet_wrap(~ period) +
    geom_point() + geom_line(aes(y=predicted))


