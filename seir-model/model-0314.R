## setwd("~/research/coronavirus/code/seir-model")

source("../configs.R")

version <- "0314"
drop.omegaeffect <- T
drop.variableomega <- F
drop.death <- T

source(paste0("modellib-", version, ".R"))

casespath <- "../../cases/panel_all.csv"
weather <- c('t2m', 'tp', 'ssrd', 'utci')
ols.priors.mu <- c(0.02773106208, -0.03134987114, -0.01027632136, -0.02558036184)
ols.priors.se <- c(0.01222135339, 0.008463339282, 0.01140204532, 0.009860867169)

df <- read.csv(casespath)
df$regid <- paste(df$Country, df$Region, df$Locality)

get.paramdf <- function(regid, param, lax, rhats, rhatparam=NULL) {
    if (is.null(lax) || length(lax) == 0)
        return(data.frame(regid, param, mu=NA, sd=NA, ci2.5=NA, ci25=NA, ci50=NA, ci75=NA, ci97.5=NA, rhat=NA))

    if (is.null(rhatparam))
        rhatparam <- param
    rhat <- mean(rhats[grep(rhatparam, rownames(rhats)), 1])
    data.frame(regid, param, mu=mean(lax), sd=sd(lax), ci2.5=quantile(lax, .025), ci25=quantile(lax, .25),
               ci50=quantile(lax, .50), ci75=quantile(lax, .75), ci97.5=quantile(lax, .975), rhat)
}

library(dplyr)
library(lfe)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

proc.model <- function(model) {
    if (drop.omegaeffect)
        model <- drop.stan.model.omega(model)
    if (drop.variableomega)
        model <- drop.stan.model.dlogomega(model)

    return(model)
}

stan.compiled <- list("full3"=list("deaths"=stan_model(model_code=proc.model(get.stan.model.deaths())),
                                   "nodice"=stan_model(model_code=proc.model(get.stan.model.nodice()))),
                      "noprior"=list("deaths"=stan_model(model_code=proc.model(drop.stan.model.prior(get.stan.model.deaths()))),
                                     "nodice"=stan_model(model_code=proc.model(drop.stan.model.prior(get.stan.model.nodice())))),
                      "noweather"=list("deaths"=stan_model(model_code=proc.model(drop.stan.model.weather(get.stan.model.deaths()))),
                                       "nodice"=stan_model(model_code=proc.model(drop.stan.model.weather(get.stan.model.nodice())))))

weatherscales <- apply(df[, weather], 2, sd)

randorder <- unique(df$regid)[sample(1:length(unique(df$regid)))]
cntyorder <- unique(df$regid[df$Region == '' & df$Locality == ''])
finalorder <- c(cntyorder, randorder[!(randorder %in% cntyorder)])

for (regid in finalorder) {
    for (model in c('noweather')) { # 'full3', 'noprior',
        subdf <- df[df$regid == regid,]

        if (regid == "Germany Berlin " && subdf$population[1] == 0)
            subdf$population <- 3769495

        outpath <- paste0("../../results/epimodel-", version, "-", model, ifelse(drop.omegaeffect, "-noomega", ""), ifelse(drop.variableomega, "-nodlogomega", ""), ifelse(drop.death, "-nodeath", ""), ".csv")

        ## Check if region is claimed
        regfile <- paste0(outpath, "-", regid)
        if (file.exists(regfile))
            next
        ## Claim this region
        fileConn <- file(regfile)
        writeLines(as.character(Sys.getpid()), fileConn)
        close(fileConn)

        print(regid)

        subdf$Confirmed[is.na(subdf$Confirmed)] <- 0
        while (sum(subdf$Confirmed[-1] < subdf$Confirmed[-nrow(subdf)]) > 0) {
            bads <- c(F, subdf$Confirmed[-1] < subdf$Confirmed[-nrow(subdf)])
            subdf$Confirmed[bads] <- c(NA, subdf$Confirmed[-nrow(subdf)])[bads]
        }

        dmobility <- diff(subdf$mobility_pca1)
        dmobility[is.na(dmobility)] <- 0 # Given the affine intercept

        stan.data <- list(T=nrow(subdf), N=round(subdf$population[1]), K=length(weather),
                          alpha_prior=.395 / 100, eein_prior=(22 / 30) / 11e6 * round(subdf$population[1]), # based on confirmed exports from Wuhan (Kucharski et al.)
                          invsigma_prior=5.2, invgamma_prior=2.9,
                          invkappa_prior=7, invtheta_prior=7,
                          beta0_prior=2.5 * 2.9, dmobility_proxy=dmobility,
                          weather=demeanlist(subdf[, weather], list(factor(rep('all', nrow(subdf))))) / t(matrix(weatherscales, ncol=nrow(subdf), nrow=length(weather))),
                          total_prior=ols.priors.mu, total_prior_sd=ols.priors.se,
                          ii_init=0, dobserved_true=diff(subdf$Confirmed) + 1)

        if (sum(!is.na(subdf$Deaths) & !is.na(subdf$Confirmed)) > 10 && !drop.death) {
            subdf$Deaths[is.na(subdf$Deaths)] <- 0
            while (sum(subdf$Deaths[-1] < subdf$Deaths[-nrow(subdf)]) > 0) {
                bads <- c(F, subdf$Deaths[-1] < subdf$Deaths[-nrow(subdf)])
                subdf$Deaths[bads] <- c(NA, subdf$Deaths[-nrow(subdf)])[bads]
            }

            stan.data$ddeaths_true <- diff(subdf$Deaths) + 1

            fit <- tryCatch({
                sampling(stan.compiled[[model]][['deaths']], data=stan.data, open_progress=F, control=list(max_treedepth=15))
            }, error=function(e) {
                NULL
            })
        } else {
            fit <- tryCatch({
                sampling(stan.compiled[[model]][['nodice']], data=stan.data, open_progress=F)
            }, error=function(e) {
                NULL
            })
        }

        la <- tryCatch({
            extract(fit, permute=T)
        }, error=function(e) {
            NULL
        })
        if (is.null(la))
            next
        rhats <- stan_rhat(fit)$data

        resrow <- rbind(get.paramdf(regid, 'alpha', la$alpha, rhats),
                        get.paramdf(regid, 'invsigma', la$invsigma, rhats),
                        get.paramdf(regid, 'invgamma', la$invgamma, rhats),
                        get.paramdf(regid, 'invkappa', la$invkappa, rhats),
                        get.paramdf(regid, 'invtheta', la$invtheta, rhats),
                        get.paramdf(regid, 'omega', la$omega, rhats),
                        get.paramdf(regid, 'mobility_slope', la$mobility_slope, rhats),
                        get.paramdf(regid, 'deathrate', la$deathrate, rhats),
                        get.paramdf(regid, 'deathlearning', la$deathlearning, rhats),
                        get.paramdf(regid, 'deathomegaplus', la$deathomegaplus, rhats),
                        get.paramdf(regid, 'error', la$error, rhats),
                        get.paramdf(regid, 'logbeta', la$logbeta, rhats),
                        get.paramdf(regid, 'logomega', la$logomega, rhats),
                        get.paramdf(regid, 'eein', la$eein / stan.data$N, rhats))
        for (kk in 1:length(weather))
            resrow <- rbind(resrow, get.paramdf(regid, paste0('e.', weather[kk]), la$effect[,kk], rhats, rhatparam=paste0('effect\\[', kk, '\\]')))
        for (kk in 1:length(weather))
            resrow <- rbind(resrow, get.paramdf(regid, paste0('o.', weather[kk]), la$omegaeffect[,kk], rhats, rhatparam=paste0('effect\\[', kk, '\\]')))

        dynrow <- data.frame(logbeta=colMeans(la$logbeta), logomega=c(NA, colMeans(la$logomega)), eein=c(NA, colMeans(la$eein)),
                             doweffect=colMeans(la$doweffect)[1 + (2:(stan.data$T+1)) %% 7],
                             dowomegaeffect=c(NA, colMeans(la$dowomegaeffect)[1 + (2:stan.data$T) %% 7]))

        write.csv(resrow, regfile, row.names=F)
        write.csv(dynrow, gsub("\\.csv", "-dynamics.csv", regfile), row.names=F)
    }
}
