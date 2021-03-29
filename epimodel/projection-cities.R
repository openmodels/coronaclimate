setwd("~/research/coronavirus/code/epimodel")

weathervars <- c('absh', 't2m', 'tp', 'ssrd', 'utci')

library(dplyr)
library(lfe)
source("forward-0105.R")
source("forward-0105-adaptive.R")

cities <- read.csv("major_cities_selection.csv")

results <- read.csv("../../results-0105/results-saved/epimodel-meta-0105noprior-all-pop.csv")

df <- read.csv("../../cases/panel_all.csv")
df$regid <- paste(df$Country, df$Region, df$Locality)

weatherscales <- apply(df[, weathervars], 2, sd)

df.global <- df %>% group_by(Country, Region, Locality) %>% summarize(Date=Date[-1], dcc=pmax(diff(Confirmed), 0), population=population[-1]) %>% group_by(Date) %>% summarize(dcc=7794798739 * sum(dcc, na.rm=T) / sum(population))

projdf <- data.frame()
logvals <- data.frame()
for (ii in 1:nrow(cities)) {
    rid <- df$regid[df$GEO_ID_REFERENCE == cities$GEO_ID_REFERENCE[ii] & df$GEO_ID == cities$GEO_ID[ii]][1]

    print(rid)
    subdf <- subset(df, regid == rid)

    for (resrid in c(rid, "  ")) {
        subres <- subset(results, regid == resrid & group == 'Combined')
        if (nrow(subres) != 24 && nrow(subres) != 23) {
            ## Country-level
            subres <- subset(results, Country == cities$Country[ii] & Region == '' & Locality == '' & group == 'Combined')
        }
        
        weather <- demeanlist(subdf[, weathervars], list(factor(rep('all', nrow(subdf))))) / t(matrix(weatherscales, ncol=nrow(subdf), nrow=length(weathervars)))

        params <- list(doweffect6=rep(0, 6), dowomegaeffect6=rep(0, 6))
        params[['logomega']] <- rep(subres$mu[subres$param == 'logomega'], nrow(weather) - 1)

        for (param in c('invsigma', 'invgamma', 'invkappa', 'invtheta',
                        'deathrate', 'deathlearning', 'deathomegaplus'))
            params[[param]] <- subres$mu[subres$param == param]
    
        params[['effect']] <- sapply(weathervars, function(var) subres$mu[subres$param == paste0('e.', var)])
        params[['omegaeffect']] <- sapply(weathervars, function(var) subres$mu[subres$param == paste0('o.', var)])

        firstcase <- which(subdf$Confirmed > 0)[1]
        params[['eein']] <- subres$mu[subres$param == 'eein'] * c(df.global$dcc[8:nrow(df.global)], rep(tail(df.global$dcc, 1), 7)) / max(df.global$dcc)
        params[['eein']][1:firstcase] <- 0
        params[['eein']][firstcase - 6] <- subdf$Confirmed[firstcase] / exp(subres$mu[subres$param == 'logomega'])

        ## Initial estimates of logbeta and logomega
        dcc <- diff(subdf$Confirmed)
        dcc[is.na(dcc)] <- 0
        fractional.shift <- c(dcc[31:length(dcc)], rep(dcc[length(dcc)], 30)) / max(dcc)
        logbeta.init <- (subres$ci97.5[subres$param == 'logbeta'] - subres$mu[subres$param == 'logbeta']) * fractional.shift + subres$mu[subres$param == 'logbeta']
        logomega.init <- (subres$ci97.5[subres$param == 'logomega'] - subres$mu[subres$param == 'logomega']) * cummax(fractional.shift) + subres$mu[subres$param == 'logomega']

        data <- list(T=nrow(weather), N=subres$population[1], K=length(weathervars), weather=weather, ii_init=0)
        
        withweather <- forward.adaptive(data, params, logbeta.init, logomega.init, subres$mu[subres$param == 'alpha'], diff(subdf$Confirmed))
        withweather$dobserved_true <- c(0, dcc)

        ## ggplot(withweather, aes(TT)) +
        ##     geom_line(aes(y=dcc, colour='Estimated')) +
        ##     geom_line(aes(y=dobserved_true, colour='Observed'))

        params[['logbeta']] <- withweather$logbeta
        params[['logomega']] <- withweather$logomega

        logvals <- rbind(logvals, data.frame(rid, logbeta=withweather$logbeta, logomega=withweather$logomega))

        ## withweather2 <- forward(data, params)

        data <- list(T=nrow(weather), N=subres$population[1], K=length(weathervars), weather=0 * weather, ii_init=0)
        baseline <- forward(data, params)
        
        subprojdf <- data.frame(regid=rid, paramregid=resrid,
                                cc0=cumsum(baseline$dcc), deaths0=cumsum(baseline$ddeaths),
                                cc1=cumsum(withweather$dcc), deaths1=cumsum(withweather$ddeaths),
                                latitude=cities$latitude[ii] * ifelse(cities$hemisphere[ii] == 'N', 1, -1))
        subprojdf$dcc <- subprojdf$cc1 - subprojdf$cc0
        subprojdf$population <- subres$population[1]
        subprojdf$maxcases <- subdf$Confirmed[nrow(subdf)]
        subprojdf$time <- 1:nrow(subprojdf)

        ## ggplot(subprojdf, aes(time, dcc / max(maxcases, cc1, na.rm=T))) +
        ##     geom_line() + geom_hline(yintercept=0, colour='green') +
        ##     theme_bw() + xlab(NULL) + ylab("Additional reported portion of total cases") + scale_y_continuous(labels=scales::percent)
        
        projdf <- rbind(projdf, subprojdf)
    }
}

library(ggplot2)

logvals.mean <- logvals %>% group_by(rid) %>% summarize(logbeta=mean(logbeta, na.rm=T), logomega=mean(logomega, na.rm=T))

projdf$date <- as.Date("2020-01-01") + projdf$time - 1

projdf.norm <- projdf %>% group_by(regid, paramregid) %>% summarize(date=date, latitude=latitude, dcc.norm=dcc / max(maxcases, cc1, na.rm=T))

ggplot(projdf.norm, aes(date, dcc.norm, linetype=paramregid == '  ')) +
    facet_grid(regid ~ ., scales='free') +
    geom_line() + geom_hline(yintercept=0, colour='green') +
    scale_x_date(expand=c(0, 0)) + scale_linetype_discrete(name="", breaks=c(F, T), labels=c('Local parameters', 'Global parameters')) +
    theme_bw() + xlab(NULL) + ylab("Additional reported portion of total cases") + scale_y_continuous(labels=scales::percent)

cities$Locality[cities$Country == 'Iceland'] <- "Reykjavík"
cities$Locality[cities$Country == 'Brazil'] <- "Brasilia"

ggplot(projdf.norm) +
    geom_hline(aes(yintercept=latitude), colour='black') +
    geom_line(aes(date, latitude + 100 * dcc.norm, group=paste(regid, paramregid), linetype=paramregid == '  ', colour=regid)) +
    geom_label(data=cities, aes(as.Date('2020-02-01'), y=latitude * ifelse(hemisphere == 'N', 1, -1), label=Locality)) +
    scale_x_date(expand=c(0, 0)) + scale_linetype_discrete(name="Calibration:", breaks=c(F, T), labels=c('Local parameters', 'Global parameters')) +
    theme_bw() + xlab(NULL) + ylab("Additional reported percent of total cases (relative to latitude)") +
    scale_colour_discrete(name="Region:")
ggsave("~/Dropbox/Coronavirus and Climate/figures/forward-cities.pdf", width=6, height=4)
