setwd("~/research/coronavirus/code/epimodel")

library(dplyr)
library(scales)
source("forward-0105.R")

results <- read.csv("../../results-saved/epimodel-meta-0105noprior-all-pop.csv")

projdf <- data.frame()
for (country in unique(results$Country)) {
    if (country == "")
        next
    subres <- subset(results, Country == country & Region == '' & Locality == '' & group == 'Combined')
    if (nrow(subres) != 24 && nrow(subres) != 23)
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
    
    for (effect in c('Transmission', 'Detection')) {
        if (effect == 'Transmission') {
            params[['effect']] <- 1
            params[['omegaeffect']] <- 0
        } else {
            params[['effect']] <- 0
            params[['omegaeffect']] <- 1
        }

        impresp <- get.impresp(baseline, 50, params)

        subprojdf <- data.frame(country, effect, time=50:(nrow(weather)-1), impresp=impresp[50:length(impresp)])

        projdf <- rbind(projdf, subprojdf)
    }
}

projdf.global <- projdf %>% group_by(effect, time) %>% summarize(country="", impresp=mean(impresp))

ggplot(subset(projdf, time <= 50 + 30), aes(time - 50, impresp)) +
    facet_grid(effect ~ .) +
    geom_line(aes(group=country), size=.2, alpha=.2) +
    geom_line(data=subset(projdf.global, time <= 50 + 30), size=1, colour=muted('red')) +
    theme_bw() + scale_x_continuous(name="Time from weather shock", expand=c(0, 0)) +
    ylab("Marginal effect of normalized weather shock")
ggsave("~/Dropbox/Coronavirus and Climate/figures/impresps.pdf", width=6, height=4)
