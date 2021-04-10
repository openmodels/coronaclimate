setwd("~/research/coronavirus/code/epimodel")

weathervars <- c('t2m', 'tp', 'ssrd', 'utci')
version <- "0314-noprior"

library(dplyr)
library(lfe)
library(ggplot2)
source("forward-0314.R")
source("forward-0314-adaptive.R")

cities <- read.csv("major_cities_selection.csv")

results <- read.csv(paste0("../../results-saved/epimodel-meta-", version, "-all-nobs-nodel.csv"))

df <- read.csv("../../cases/panel_all.csv")
df$regid <- paste(df$Country, df$Region, df$Locality)

weatherscales <- apply(df[, weathervars], 2, sd)

df.global <- df %>% group_by(Country, Region, Locality) %>% summarize(Date=Date[-1], dcc=pmax(diff(Confirmed), 0), population=population[-1]) %>% group_by(Date) %>% summarize(dcc=7794798739 * sum(dcc, na.rm=T) / sum(population))

projdf <- data.frame()
sumstats <- data.frame()
for (ii in 1:nrow(cities)) {
    ridinfo <- df[df$GEO_ID_REFERENCE == cities$GEO_ID_REFERENCE[ii] & df$GEO_ID == cities$GEO_ID[ii],][1,]
    print(ridinfo$regid)

    ridlevels <- data.frame(Country=c(rep(ridinfo$Country, 3), ""), Region=c(rep(ridinfo$Region, 2), rep("", 2)), Locality=c(ridinfo$Locality, rep("", 3)))
    ridlevels <- ridlevels[!duplicated(ridlevels),]
    rids <- paste(ridlevels$Country, ridlevels$Region, ridlevels$Locality)
    
    subdf <- subset(df, regid == ridinfo$regid)

    for (global.params in c(F, T)) {
        if (global.params) {
            paramregid <- "  "
            subres <- subset(results, regid == paramregid & group == 'Combined')
        } else {
            for (paramregid in rids) {
                subres <- subset(results, regid == paramregid & group == 'Combined')
                if (nrow(subres) == 22 || nrow(subres) == 21)
                    break
            }
        }
        
        weather <- demeanlist(subdf[, weathervars], list(factor(rep('all', nrow(subdf))))) / t(matrix(weatherscales, ncol=nrow(subdf), nrow=length(weathervars)))

        for (dynregid in rids) {
            dynamic.paths <- Sys.glob(paste0("../../results-saved/epimodel-", version, "-dynamics.csv-", gsub(" +$", "*", dynregid)))
            if (length(dynamic.paths) > 0) {
                dynamics <- NULL
                for (dynamic.path in dynamic.paths) {
                    dynamics.one <- read.csv(dynamic.path)
                    if (is.null(dynamics))
                        dynamics <- dynamics.one
                    else
                        dynamics <- dynamics + dynamics.one
                }
                dynamics <- dynamics / length(dynamic.paths)
                break
            }
        }
        if (dynregid == "  ") {
            print("No dynamics found.")
            next
        }
        
        params <- list(doweffect6=c(dynamics$doweffect[6:7], dynamics$doweffect[1:4]),
                       dowomegaeffect6=c(dynamics$dowomegaeffect[7:8], dynamics$dowomegaeffect[2:5]),
                       logbeta=dynamics$logbeta, logomega=dynamics$logomega[-1], eein=dynamics$eein[-1])

        for (param in c('invsigma', 'invgamma', 'invkappa', 'invtheta',
                        'deathrate', 'deathlearning', 'deathomegaplus'))
            params[[param]] <- subres$mu[subres$param == param]
    
        params[['effect']] <- sapply(weathervars, function(var) subres$mu[subres$param == paste0('e.', var)])
        params[['omegaeffect']] <- sapply(weathervars, function(var) subres$mu[subres$param == paste0('o.', var)])

        data <- list(T=nrow(weather), N=subres$population[1], K=length(weathervars), weather=weather, ii_init=0)
        
        withweather <- forward.adaptive(data, params, diff(subdf$Confirmed))
        withweather$dobserved_true <- c(0, diff(subdf$Confirmed))

        ## ggplot(withweather, aes(TT)) +
        ##     geom_line(aes(y=dcc, colour='Estimated')) +
        ##     geom_line(aes(y=dobserved_true, colour='Observed'))

        rsqr <- 1 - sum((withweather$dobserved_true - withweather$dcc)^2, na.rm=T) / sum(withweather$dobserved_true^2, na.rm=T)
        sumstats <- rbind(sumstats, data.frame(regid=ridinfo$regid, dynregid, paramregid, rsqr))
        
        ## withweather2 <- forward(data, params, withweather$extraeein)
        ## withweather2$dobserved_true <- c(0, diff(subdf$Confirmed))
        ## ggplot(withweather2, aes(TT)) +
        ##     geom_line(aes(y=dcc, colour='Estimated')) +
        ##     geom_line(aes(y=dobserved_true, colour='Observed'))

        data <- list(T=nrow(weather), N=subres$population[1], K=length(weathervars), weather=0 * weather, ii_init=0)
        baseline <- forward(data, params, withweather$extraeein)
        
        subprojdf <- data.frame(regid=ridinfo$regid, dynregid, paramregid, global.params,
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

projdf$date <- as.Date("2020-01-01") + projdf$time - 1

projdf.norm <- projdf %>% group_by(regid, paramregid) %>% summarize(global.params, date=date, latitude=latitude, dcc.norm=dcc / max(maxcases, cc0, cc1, na.rm=T))

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

ggplot(subset(projdf.norm, !global.params)) +
    geom_hline(aes(yintercept=latitude, colour=regid), linetype='dashed') +
    geom_line(aes(date, latitude + 100 * dcc.norm, group=paste(regid, paramregid), colour=regid)) +
    geom_label(data=cities, aes(as.Date('2020-02-01'), y=latitude * ifelse(hemisphere == 'N', 1, -1), label=Locality)) +
    scale_x_date(expand=c(0, 0)) + scale_linetype_discrete(name="Calibration:", breaks=c(F, T), labels=c('Local parameters', 'Global parameters')) +
    theme_bw() + xlab(NULL) + ylab("Additional reported percent of total cases (relative to latitude)") +
    scale_colour_discrete(name="Region:")

ggsave("~/Dropbox/Coronavirus and Climate/figures/forward-cities-0314.pdf", width=8, height=6)
