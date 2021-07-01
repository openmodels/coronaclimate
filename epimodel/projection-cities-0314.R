setwd("~/research/coronavirus/code/epimodel")

weathervars <- c('t2m', 'tp', 'ssrd', 'utci')
version <- "0314-noprior"

library(dplyr)
library(lfe)
library(ggplot2)
source("forward-0314.R")
source("forward-0314-adaptive.R")

cities <- read.csv("major_cities_selection.csv")

results <- read.csv(paste0("../../results/epimodel-meta-", version, "-all-nobs-nodel.csv"))

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
            dynamic.paths <- Sys.glob(paste0("../../results/epimodel-", version, "-dynamics.csv-", gsub(" +$", "*", dynregid)))
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

        data <- list(T=nrow(weather), N=subdf$population[1], K=length(weathervars), weather=weather, ii_init=0)

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

        print(c(ridinfo$regid, dynregid, paramregid))
        
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

projdf.norm <- projdf %>% group_by(regid, paramregid) %>% summarize(global.params, date=date, latitude=latitude, dcc.norm=pmax(dcc / max(maxcases, cc0, cc1, na.rm=T), -1))

gp <- ggplot(projdf.norm[projdf.norm$global.params == F,], aes(date, dcc.norm, linetype=paramregid == '  ')) +
    facet_grid(regid ~ ., scales='free') +
    geom_line() + geom_hline(yintercept=0, colour='green') +
    scale_x_date(expand=c(0, 0)) + scale_linetype_discrete(name="", breaks=c(F, T), labels=c('Local parameters', 'Global parameters')) +
    theme_bw() + xlab(NULL) + ylab("Additional reported portion of total cases") + scale_y_continuous(labels=scales::percent)
ggsave("tmp.pdf", gp, width=8, height=8)

cities$Locality[cities$Country == 'Iceland'] <- "Reykjavík"
cities$Locality[cities$Country == 'Brazil'] <- "Sao Paulo"

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

## Add in internal axes
verticals <- data.frame()
ticks <- data.frame()
for (ii in 1:nrow(cities)) {
    latitude <- cities$latitude[ii] * ifelse(cities$hemisphere[ii] == 'N', 1, -1)
    rows <- which(projdf.norm$latitude == latitude & !projdf.norm$global.params)
    ymin <- latitude + 100 * min(projdf.norm$dcc.norm[rows], na.rm=T)
    ymax <- latitude + 100 * max(projdf.norm$dcc.norm[rows], na.rm=T)

    ##date <- max(max(projdf.norm$date) - 10, verticals$date[verticals$ymin < ymax + 10 & verticals$ymax > ymin - 10]) + 22
    for (date in as.character(seq(max(projdf.norm$date) + 12, by=22, length.out=10))) {
        date <- as.Date(date)
        if (!any(verticals$date == date & verticals$ymin < ymax + 10 & verticals$ymax > ymin - 10))
            break
    }
    verticals <- rbind(verticals, data.frame(date, latitude, ymin, ymax, topval=ymax - latitude, botval=ymin - latitude))

    values <- c(seq(0, 100 * max(projdf.norm$dcc.norm[rows], na.rm=T), by=10),
                seq(0, 100 * min(projdf.norm$dcc.norm[rows], na.rm=T), by=-10)[-1])
    ticks <- rbind(ticks, data.frame(y=latitude + values, date, values, latitude))
}

gp <- ggplot(subset(projdf.norm, !global.params)) +
    geom_segment(data=verticals, aes(x=as.Date('2020-01-01'), xend=date, y=latitude, yend=latitude, colour=factor(latitude)), linetype='dashed') +
    geom_line(aes(date, latitude + 100 * dcc.norm, group=paste(regid, paramregid), colour=factor(latitude))) +
    geom_label(data=cities, aes(as.Date('2020-02-01'), y=latitude * ifelse(hemisphere == 'N', 1, -1), label=Locality)) +
    geom_segment(data=verticals, aes(x=date, xend=date, y=ymin, yend=ymax, colour=factor(latitude))) +
    geom_segment(data=ticks, aes(x=date, xend=date+4, y=y, yend=y, colour=factor(latitude))) +
    geom_text(data=subset(ticks, values == 0), aes(x=date+6, y=y, label=values), hjust=0) +
    geom_text(data=subset(verticals, round(botval) != 0), aes(x=date, y=ymin - 2, label=paste0(round(botval), "%"))) +
    geom_text(data=subset(verticals, round(topval) != 0), aes(x=date, y=ymax + 2, label=paste0(round(topval), "%"))) +
    scale_x_date(expand=c(0, 0)) + scale_linetype_discrete(name="Calibration:", breaks=c(F, T), labels=c('Local parameters', 'Global parameters')) +
    theme_bw() + xlab(NULL) + ylab("Additional reported percent of total cases (relative to latitude)") +
    guides(colour=F) + coord_cartesian(xlim=c(as.Date("2020-01-01"), as.Date("2020-12-31")), clip="off") + theme(plot.margin=margin(5, 150, 5, 5, "pt"))

ggsave("~/Dropbox/Coronavirus and Climate/figures/forward-cities-0314.pdf", gp, width=9, height=6)

## Do all regions

sumstats <- data.frame()
for (regid in unique(results$regid[results$lowest_level == 1 & results$group == "Combined"])) {
    print(regid)
    subdf <- df[df$regid == regid,]
    subres <- results[results$regid == regid & results$group == 'Combined',]
    if (nrow(subres) != 22 && nrow(subres) != 21)
        next

    weather <- demeanlist(subdf[, weathervars], list(factor(rep('all', nrow(subdf))))) / t(matrix(weatherscales, ncol=nrow(subdf), nrow=length(weathervars)))

    dynamics <- read.csv(paste0("../../results/epimodel-", version, "-dynamics.csv-", regid))

    params <- list(doweffect6=c(dynamics$doweffect[6:7], dynamics$doweffect[1:4]),
                   dowomegaeffect6=c(dynamics$dowomegaeffect[7:8], dynamics$dowomegaeffect[2:5]),
                   logbeta=dynamics$logbeta, logomega=dynamics$logomega[-1], eein=dynamics$eein[-1])

    for (param in c('invsigma', 'invgamma', 'invkappa', 'invtheta',
                    'deathrate', 'deathlearning', 'deathomegaplus'))
        params[[param]] <- subres$mu[subres$param == param]

    params[['effect']] <- sapply(weathervars, function(var) subres$mu[subres$param == paste0('e.', var)])
    params[['omegaeffect']] <- sapply(weathervars, function(var) subres$mu[subres$param == paste0('o.', var)])

    data <- list(T=nrow(weather), N=subres$population[1], K=length(weathervars), weather=weather, ii_init=0)

    withweather <- tryCatch({
        forward.adaptive(data, params, diff(subdf$Confirmed))
    }, error=function(e) {
        NULL
    })
    if (is.null(withweather))
        next
    withweather$dobserved_true <- c(0, diff(subdf$Confirmed))

    rsqr <- 1 - sum((withweather$dobserved_true - withweather$dcc)^2, na.rm=T) / sum(withweather$dobserved_true^2, na.rm=T)

    data <- list(T=nrow(weather), N=subres$population[1], K=length(weathervars), weather=0 * weather, ii_init=0)
    baseline <- forward(data, params, withweather$extraeein)

    subprojdf <- data.frame(regid,
                            cc0=cumsum(baseline$dcc), deaths0=cumsum(baseline$ddeaths),
                            cc1=cumsum(withweather$dcc), deaths1=cumsum(withweather$ddeaths),
                            latitude=cities$latitude[ii] * ifelse(cities$hemisphere[ii] == 'N', 1, -1))
    subprojdf$dcc <- subprojdf$cc1 - subprojdf$cc0
    subprojdf$ddeaths <- subprojdf$deaths1 - subprojdf$deaths0
    subprojdf$time <- 1:nrow(subprojdf)

    population <- subres$population[1]
    maxcases <- subdf$Confirmed[nrow(subdf)]
    maxdeaths <- subdf$Deaths[nrow(subdf)]

    subprojdf$dcc.norm <- subprojdf$dcc / max(maxcases, subprojdf$cc0, subprojdf$cc1, na.rm=T)
    subprojdf$ddeaths.norm <- subprojdf$ddeaths / max(maxdeaths, subprojdf$deaths0, subprojdf$deaths1, na.rm=T)

    sumstats <- rbind(sumstats, data.frame(regid, rsqr, population, maxcases, max.dcc.norm=max(subprojdf$dcc.norm, na.rm=T), min.dcc.norm=min(subprojdf$dcc.norm, na.rm=T), end.dcc.norm=tail(subprojdf$dcc.norm, 1), max.ddeaths.norm=max(subprojdf$ddeaths.norm, na.rm=T), min.ddeaths.norm=min(subprojdf$ddeaths.norm, na.rm=T), end.ddeaths.norm=tail(subprojdf$ddeaths.norm, 1), var.dcc0=var(baseline$dcc / max(maxcases, subprojdf$cc0, subprojdf$cc1, na.rm=T), na.rm=T), var.dcc1=var(withweather$dcc / max(maxcases, subprojdf$cc0, subprojdf$cc1, na.rm=T), na.rm=T)))
}

library(Hmisc)
mean(sumstats$var.dcc1) / mean(sumstats$var.dcc0)
wtd.mean(sumstats$var.dcc1, sumstats$population) / wtd.mean(sumstats$var.dcc0, sumstats$population)

ggplot(sumstats) +
    geom_histogram(aes(max.dcc.norm, y=..density.., weight=population, fill='Maximum')) + geom_histogram(aes(min.dcc.norm, y=..density.., weight=population, fill='Minimum'))

sumstats$ext.dcc.norm <- sumstats$max.dcc.norm
sumstats$ext.dcc.norm[-sumstats$min.dcc.norm > sumstats$ext.dcc.norm] <- sumstats$min.dcc.norm[-sumstats$min.dcc.norm > sumstats$ext.dcc.norm]

sumstats$ext.ddeaths.norm <- sumstats$max.ddeaths.norm
sumstats$ext.ddeaths.norm[-sumstats$min.ddeaths.norm > sumstats$ext.ddeaths.norm] <- sumstats$min.ddeaths.norm[-sumstats$min.ddeaths.norm > sumstats$ext.ddeaths.norm]

sumstats$bin <- round(2 * sumstats$ext.dcc.norm, 1) / 2
sumstats.ext.dcc <- sumstats %>% group_by(bin) %>% dplyr::summarize(stat='Peak Cases', dpop=sum(population))
sumstats$bin <- round(2 * sumstats$end.dcc.norm, 1) / 2
sumstats.end.dcc <- sumstats %>% group_by(bin) %>% dplyr::summarize(stat='Year-end Cases', dpop=sum(population))
sumstats$bin <- round(2 * sumstats$ext.ddeaths.norm, 1) / 2
sumstats.ext.ddeaths <- sumstats %>% group_by(bin) %>% dplyr::summarize(stat='Peak Deaths', dpop=sum(population))
sumstats$bin <- round(2 * sumstats$end.ddeaths.norm, 1) / 2
sumstats.end.ddeaths <- sumstats %>% group_by(bin) %>% dplyr::summarize(stat='Year-end Deaths', dpop=sum(population))

sumstats.agg <- rbind(sumstats.ext.dcc, sumstats.end.dcc, sumstats.ext.ddeaths, sumstats.end.ddeaths)

gp <- ggplot(sumstats.agg) +
    facet_wrap(~ stat, ncol=2, nrow=2, scales="free_y") +
    geom_col(aes(-bin, dpop)) + # flip sign, so it's cc0 - cc1, consistent with normalization
    scale_x_continuous(expand=c(0, 0), labels=scales::percent) + scale_y_continuous(expand=expansion(mult=c(0, .05))) +
    theme_bw() + xlab("Change in the absence of weather variation") +
    ylab("Populations of regions affected")
ggsave("~/Dropbox/Coronavirus and Climate/figures/forward-cities-0314-hists.pdf", gp, width=9, height=6)

gp <- ggplot(subset(sumstats.agg, stat %in% c('Peak Cases', 'Peak Deaths'))) +
    facet_wrap(~ stat, ncol=2, scales="free_y") +
    geom_col(aes(-bin, dpop)) + # flip sign, so it's cc0 - cc1, consistent with normalization
    scale_x_continuous(expand=c(0, 0), labels=scales::percent) + scale_y_continuous(expand=expansion(mult=c(0, .05))) +
    theme_bw() + xlab("Change in the absence of weather variation") +
    ylab("Populations of regions affected")
ggsave("~/Dropbox/Coronavirus and Climate/figures/forward-cities-0314-peakhists.pdf", gp, width=9, height=3)

sumstats$bin <- round(4 * sumstats$max.dcc.norm, 1) / 4
sumstats.max.dcc <- sumstats %>% group_by(bin) %>% dplyr::summarize(stat='Maximum Cases Change', dpop=sum(population))
sumstats$bin <- round(4 * sumstats$max.ddeaths.norm, 1) / 4
sumstats.max.ddeaths <- sumstats %>% group_by(bin) %>% dplyr::summarize(stat='Maximum Deaths Change', dpop=sum(population))

sumstats.agg <- rbind(sumstats.max.dcc, sumstats.max.ddeaths)

gp <- ggplot(sumstats.agg) +
    facet_wrap(~ stat, nrow=1, scales="free_y") +
    geom_col(aes(bin, dpop)) +
    scale_x_continuous(expand=c(0, 0), labels=scales::percent) + scale_y_continuous(expand=expansion(mult=c(0, .05))) +
    theme_bw() + xlab("Change relative to constant weather") +
    ylab("Populations of regions affected")

ggsave("~/Dropbox/Coronavirus and Climate/figures/forward-cities-0314-maxhists.pdf", gp, width=9, height=3)

sum(sumstats$population[sumstats$max.dcc.norm > .21]) / sum(sumstats$population)
