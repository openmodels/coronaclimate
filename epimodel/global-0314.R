## setwd("~/research/coronavirus/code/epimodel")

source("../configs.R")

version <- "0314"

source(paste0("modellib-", version, ".R"))

casespath <- "../../cases/panel_all_fixed.csv"

df <- read.csv(casespath)
df$regid <- paste(df$Country, df$Region, df$Locality)

library(dplyr)

subdf <- subset(df, lowest_level == 1) %>% group_by(Date) %>%
    summarize(Country="", Region="", Locality="", Confirmed=sum(Confirmed, na.rm=T) / 1e6, Deaths=sum(Deaths, na.rm=T) / 1e6,
              Recovered=sum(Recovered, na.rm=T) / 1e6, Source="multiple", population=sum(population, na.rm=T) / 1e6)

subdf$Confirmed[is.na(subdf$Confirmed)] <- 0
while (sum(subdf$Confirmed[-1] < subdf$Confirmed[-nrow(subdf)]) > 0) {
    bads <- c(F, subdf$Confirmed[-1] < subdf$Confirmed[-nrow(subdf)])
    subdf$Confirmed[bads] <- c(NA, subdf$Confirmed[-nrow(subdf)])[bads]
}

dmobility <- rep(0, nrow(subdf)-1)
print(subdf$population[1])

stan.data <- list(T=nrow(subdf), N=round(subdf$population[1]),
                  alpha_prior=.395 / 100, eein_prior=(22 / 30) / 11e6 * round(subdf$population[1]),
                  invsigma_prior=5.2, invgamma_prior=2.9,
                  invkappa_prior=7, invtheta_prior=7,
                  beta0_prior=2.5 * 2.9, dmobility_proxy=dmobility,
                  ii_init=0, dobserved_true=diff(subdf$Confirmed) + 1)

if (!file.exists("global-0314.RData")) {
    library(rstan)
    options(mc.cores = parallel::detectCores())
    rstan_options(auto_write = TRUE)

    stan.noweather.nodice <- stan_model(model_code=drop.stan.model.weather(get.stan.model.nodice()))

    fit <- sampling(stan.noweather.nodice, data=stan.data, open_progress=F)
    print(fit)

    la <- extract(fit, permute=T)
    rhats <- stan_rhat(fit)$data

    get.paramdf <- function(param, lax, rhats, rhatparam=NULL) {
        if (is.null(lax) || length(lax) == 0)
            return(data.frame(param, mu=NA, sd=NA, ci2.5=NA, ci25=NA, ci50=NA, ci75=NA, ci97.5=NA, rhat=NA))
        
        if (is.null(rhatparam))
            rhatparam <- param
        rhat <- mean(rhats[grep(rhatparam, rownames(rhats)), 1])
        data.frame(param, mu=mean(lax), sd=sd(lax), ci2.5=quantile(lax, .025), ci25=quantile(lax, .25),
                   ci50=quantile(lax, .50), ci75=quantile(lax, .75), ci97.5=quantile(lax, .975), rhat)
    }

    resrow <- rbind(get.paramdf('alpha', la$alpha, rhats),
                    get.paramdf('invsigma', la$invsigma, rhats),
                    get.paramdf('invgamma', la$invgamma, rhats),
                    get.paramdf('invkappa', la$invkappa, rhats),
                    get.paramdf('invtheta', la$invtheta, rhats),
                    get.paramdf('omega', la$omega, rhats),
                    get.paramdf('mobility_slope', la$mobility_slope, rhats),
                    get.paramdf('deathrate', la$deathrate, rhats),
                    get.paramdf('deathlearning', la$deathlearning, rhats),
                    get.paramdf('deathomegaplus', la$deathomegaplus, rhats),
                    get.paramdf('error', la$error, rhats),
                    get.paramdf('logbeta', la$logbeta, rhats),
                    get.paramdf('logomega', la$logomega, rhats),
                    get.paramdf('eein', la$eein / stan.data$N, rhats))

    print(resrow)

    save(la, rhats, file="global-0314.RData")
}

load("global-0314.RData")

weatherscales <- apply(df[, weather], 2, sd)

source("forward-0314.R")

stan.data$K <- 1
stan.data$weather <- matrix(rep(0, nrow(subdf)), nrow(subdf), 1)

results <- read.csv(paste0("../../results/epimodel-meta-", version, "-noprior-all-nobs-nodel.csv"))

subres <- subset(results, regid == "  " & group == 'Combined')
params.all <- list(logbeta=colMeans(la$logbeta), logomega=colMeans(la$logomega), eein=colMeans(la$eein))

for (param in c('invsigma', 'invgamma', 'invkappa', 'invtheta',
                'deathrate', 'deathlearning', 'deathomegaplus'))
    params.all[[param]] <- subres$mu[subres$param == param]

params.all[['effect']] <- subres$mu[subres$param == paste0('e.t2m')]
params.all[['omegaeffect']] <- subres$mu[subres$param == paste0('o.t2m')]

primed.all <- list(doweffect=colMeans(la$doweffect), dowomegaeffect=colMeans(la$dowomegaeffect),
                   ss=colMeans(la$ss), ee1=colMeans(la$ee1), ee2=colMeans(la$ee2), ii1=colMeans(la$ii1),
                   ii2=colMeans(la$ii2), qq=colMeans(la$qq), rr=colMeans(la$rr))

allprojs <- data.frame()
for (tt in 1:364) {
    print(tt)
    data <- stan.data
    data$T <- 366 - tt

    primed <- primed.all
    for (name in names(primed)) {
        if (name %in% c('doweffect', 'dowomegaeffect'))
            primed[[name]] <- c(primed[[name]][-1], primed[[name]][1])
        else
            primed[[name]] <- primed[[name]][tt:length(primed[[name]])]
    }

    params <- params.all
    params$logbeta <- params$logbeta[tt:length(params$logbeta)]
    params$logomega <- params$logomega[tt:length(params$logomega)]
    params$eein <- params$eein[tt:length(params$eein)]
    
    projs0 <- forward(data, params, primed=primed)

    data$weather <- matrix(c(rep(1, 15), rep(0, nrow(subdf) - 15)) / weatherscales[weather == 't2m'], nrow(subdf), 1)
    
    projs1 <- forward(data, params, primed=primed)

    allprojs <- rbind(allprojs, data.frame(Date=subdf$Date[tt:nrow(subdf)], diff=c(0, cumsum(projs1$dcc - projs0$dcc)),
                                           C0pdiff=subdf$Confirmed[tt] + c(0, cumsum(projs1$dcc - projs0$dcc)), start=tt,
                                           Confirmed=subdf$Confirmed[tt:nrow(subdf)] + c(0, cumsum(projs1$dcc - projs0$dcc))))
}

library(ggplot2)

gp <- ggplot(allprojs, aes(as.Date(Date), C0pdiff, group=factor(start), colour=start)) +
    geom_line() + scale_x_date() + xlab(NULL) + ylab("Confirmed")
ggsave("~/Dropbox/Coronavirus and Climate/figures/global-project-0314.pdf", gp, width=5, height=5)

range(allprojs$diff)

gp <- ggplot(subdf, aes(as.Date(Date), Confirmed)) +
    geom_line() + geom_line(data=subset(allprojs, start == allprojs$start[allprojs$diff == min(allprojs$diff)]), col='red') +
    scale_x_date()
ggsave("~/Dropbox/Coronavirus and Climate/figures/global-project-0314-max.pdf", gp, width=5, height=5)
