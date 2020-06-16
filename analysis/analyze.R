setwd("~/Dropbox/Coronavirus and Climate")

source("code/analysis/load.R")

mod <- felm(dlog ~ t2m + t2m2 + tp | factor(regid) + factor(regid) : days + factor(superset) : factor(Date) | 0 | regid, data=df)
summary(mod)

mod <- felm(dlog ~ t2m + t2m2 + tp + tp2 | factor(regid) + factor(regid) : factor(week) + factor(superset) : factor(Date) | 0 | regid, data=df)
summary(mod)

mod <- felm(dlog ~ absh + ssrd + t2m + t2m2 + tp | factor(regid) + factor(regid) : factor(week) + factor(superset) : factor(Date) | 0 | regid, data=df)
summary(mod)

source("code/analysis/lib.R")

preddf <- data.frame(tmean=seq(-10, 40), tp=0, tp2=0)
preddf$t2m <- preddf$tmean - 20
preddf$t2m2 <- preddf$tmean^2 - 20^2
preddf$t2m3 <- preddf$tmean^3 - 20^3

plot.doseresp(preddf, "Average temperature", "Change in growth rate", indep='tmean')

ggplot(df, aes(tas)) +
    geom_histogram()

## One country at a time
mods <- list("global"=mod)

for (country in names(table(df$Country)[table(df$Country) > 85])) {
    print(country)
    mod <- felm(dlog ~ tas + tas2 + tas3 + tp | factor(regid) + factor(regid) : days + factor(Date) | 0 | regid, data=subset(df, Country == country))
    summary(mod)

    mods[[country]] <- mod
}

allplotdf <- data.frame()
for (country in names(mods)) {
    plotdf <- predict.felm(mods[[country]], preddf, interval="confidence")
    plotdf$country <- country
    plotdf$tmean <- preddf$tmean

    allplotdf <- rbind(allplotdf, plotdf)
}

ggplot(allplotdf[allplotdf$country != 'global',], aes(tmean, fit)) +
    facet_wrap(~ country) +
    geom_line() + geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=.5) + theme_bw() +
    scale_x_continuous(expand=c(0, 0)) + xlab("Average temperature") + ylab("Change in growth rate") +
    scale_y_continuous(labels=scales::percent) +
    coord_cartesian(ylim=c(-2, 2))

