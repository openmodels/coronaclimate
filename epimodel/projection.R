setwd("~/research/coronavirus/code/epimodel")

weathervars <- c('t2m', 'tp', 'ssrd', 'utci')

library(lfe)
source("forward-0105.R")
source("forward-0105-adaptive.R")

results <- read.csv("../../results-saved/epimodel-meta-0105noprior-all-pop.csv")

df <- read.csv("../../cases/panel_all.csv")
weatherscales <- apply(df[, weathervars], 2, sd)

projdf <- data.frame()
for (country in unique(df$Country)) {
    subres <- subset(results, Country == country & Region == '' & Locality == '' & group == 'Combined')
    if (nrow(subres) != 24)
        next
    
    print(country)
    subdf <- subset(df, Country == country & Region == '' & Locality == '')
    
    weather <- demeanlist(subdf[, weathervars], list(factor(rep('all', nrow(subdf))))) / t(matrix(weatherscales, ncol=nrow(subdf), nrow=length(weathervars)))

    params <- list(doweffect6=rep(0, 6), dowomegaeffect6=rep(0, 6))
    params[['logomega']] <- rep(subres$mu[subres$param == 'logomega'], nrow(weather) - 1)

    for (param in c('invsigma', 'invgamma', 'invkappa', 'invtheta',
                    'deathrate', 'deathlearning', 'deathomegaplus'))
        params[[param]] <- subres$mu[subres$param == param]
    
    params[['effect']] <- sapply(weathervars, function(var) subres$mu[subres$param == paste0('e.', var)])
    params[['omegaeffect']] <- sapply(weathervars, function(var) subres$mu[subres$param == paste0('o.', var)])

    firstcase <- which(subdf$Confirmed > 0)[1]
    params[['eein']] <- rep(0, nrow(weather) - 1) # subres$mu[subres$param == 'eein']
    params[['eein']][firstcase - 6] <- subdf$Confirmed[firstcase]

    data <- list(T=nrow(weather), N=subres$population[1], K=length(weathervars), weather=weather, ii_init=0)
    withweather <- forward.adaptive(data, params, subres$mu[subres$param == 'logbeta'], subres$mu[subres$param == 'logomega'], subres$mu[subres$param == 'alpha'], diff(subdf$Confirmed))

    params[['logbeta']] <- withweather$logbeta
    params[['logomega']] <- withweather$logomega

    ## withweather2 <- forward(data, params)

    data <- list(T=nrow(weather), N=subres$population[1], K=length(weathervars), weather=0 * weather, ii_init=0)
    baseline <- forward(data, params)

    subprojdf <- data.frame(country, cc0=cumsum(baseline$dcc), deaths0=cumsum(baseline$ddeaths),
                            cc1=cumsum(withweather$dcc), deaths1=cumsum(withweather$ddeaths))
    subprojdf$dcc <- subprojdf$cc1 - subprojdf$cc0
    subprojdf$population <- subres$population[1]
    subprojdf$maxcases <- subdf$Confirmed[nrow(subdf)]
    subprojdf$time <- 1:nrow(subprojdf)

    projdf <- rbind(projdf, subprojdf)
}

plot(subdf$Confirmed, col=3)
lines(subprojdf$cc0)
lines(subprojdf$cc1, col=2)

plot(diff(subdf$Confirmed), col=3)
points(baseline$dcc)

plot(baseline$logomega)
lines(baseline$logbeta, col=2)

library(ggplot2)

projdf$date <- as.Date("2020-01-01") + projdf$time - 1

ggplot(projdf, aes(date, dcc / max(maxcases, cc1, na.rm=T))) +
    geom_line(aes(group=country), alpha=.33) + geom_hline(yintercept=0, colour='green') + geom_smooth() +
    scale_x_date(expand=c(0, 0)) +
    theme_bw() + xlab(NULL) + ylab("Additional reported portion of total cases") + scale_y_continuous(labels=scales::percent)
