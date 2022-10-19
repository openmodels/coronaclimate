setwd("~/research/coronavirus/code/epimodel")

library(dplyr)
library(scales)
source("forward-0314.R")

results <- read.csv("../../results-saved/epimodel-meta-0314-full-all-nobs-nodel.csv")

projdf <- data.frame()
for (country in unique(results$Country)) {
    if (country == "")
        next
    subres <- subset(results, Country == country & Region == '' & Locality == '' & group == 'Combined')
    if (nrow(subres) < 15 || nrow(subres) > 24)
        next
    
    print(country)

    weather <- matrix(0, 365, 1)
    data <- list(T=365, N=1, K=1, weather=weather, ii_init=1 / 10e6)

    params <- list(doweffect6=rep(0, 6), dowomegaeffect6=rep(0, 6))
    params[['logbeta']] <- rep(subres$mu[subres$param == 'logbeta'], nrow(weather))
    params[['logomega']] <- rep(subres$mu[subres$param == 'logomega'], nrow(weather) - 1)
    params[['eein']] <- rep(0, nrow(weather) - 1)

    for (param in c('invsigma', 'invgamma', 'invkappa', 'invtheta',
                    'deathrate', 'deathlearning', 'deathomegaplus'))
        params[[param]] <- subres$mu[subres$param == param]

    params[['effect']] <- 0
    params[['omegaeffect']] <- 0

    baseline <- get.dlog(data, params)

    baseline[!is.finite(baseline)] <- NA
    maxgrowth <- which.max(baseline)
    
    for (effect in c('Transmission', 'Detection')) {
        if (effect == 'Transmission') {
            params[['effect']] <- 1
            params[['omegaeffect']] <- 0
        } else {
            params[['effect']] <- 0
            params[['omegaeffect']] <- 1
        }

        impresp <- get.impresp(baseline, maxgrowth, params)

        subprojdf <- data.frame(country, effect, time=0:(nrow(weather)-maxgrowth-1), impresp=impresp[maxgrowth:length(impresp)])

        projdf <- rbind(projdf, subprojdf)
    }
}

projdf.global <- projdf %>% group_by(effect, time) %>% summarize(country="", impresp=mean(impresp))

ggplot(subset(projdf, time <= 30)) +
    facet_grid(effect ~ .) + coord_cartesian(ylim=c(-.05, .1)) +
    geom_line(aes(time, impresp, group=country), size=.2, alpha=.2) +
    geom_line(data=subset(projdf.global, time <= 30), aes(time, impresp), size=1, colour=muted('red')) +
    geom_rect(data=data.frame(effect=c('Transmission', 'Detection')), aes(xmin=0, xmax=3.5, ymin=-.1, ymax=.2, fill='one'), alpha=.25) +
    geom_rect(data=data.frame(effect=c('Transmission', 'Detection')), aes(xmin=3.5, xmax=15.5, ymin=-.1, ymax=.2, fill='two'), alpha=.25) +
    geom_rect(data=data.frame(effect=c('Transmission', 'Detection')), aes(xmin=15.5, xmax=30, ymin=-.1, ymax=.2, fill='three'), alpha=.25) +
    theme_bw() + scale_x_continuous(name="Time from weather shock", expand=c(0, 0)) +
    ylab("Marginal effect of normalized weather shock") +
    scale_fill_discrete(name="Periods: ", breaks=c('one', 'two', 'three'), labels=c('Detection-only', 'Counter effects', 'Transmission-only')) +
    theme(legend.position="bottom")

ggsave("~/Dropbox/Coronavirus and Climate/figures/impresps.pdf", width=6, height=4.5)
