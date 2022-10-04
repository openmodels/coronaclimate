## setwd("~/research/coronavirus/code/epimodel")

source("../configs.R")

source("forward-0314.R")

## Generate simulated data

weather <- matrix(-cos(2*pi*(0:364)/365), 365, 1)

data <- list(T=365, N=1e6, K=1, weather=weather, ii_init=1)

params <- list(invsigma=5.2, invgamma=2.9, invkappa=7, invtheta=7,
               beta0=2.5 * 2.9, doweffect6=rep(0, 6), dowomegaeffect6=rep(0, 6),
               logbeta=0.5 * ((0:364) / 365) - 0.7, logomega=0.5 * ((0:364) / 365) - 0.7,
               effect=-.01, omegaeffect=.01, deathrate=.01, deathlearning=.01, deathomegaplus=.2,
               eein=rep(0, 365))

out <- forward(data, params)

plot(out$TT, out$dcc)

version <- "0314"

source(paste0("modellib-", version, ".R"))

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Generate simulated data

stan.data <- list(T=365, N=1e6, K=1,
                  alpha_prior=.395 / 100, eein_prior=(22 / 30) / 11e6 * 1e6, # based on confirmed exports from Wuhan (Kucharski et al.)
                  invsigma_prior=5.2, invgamma_prior=2.9,
                  invkappa_prior=7, invtheta_prior=7,
                  beta0_prior=2.5 * 2.9, dmobility_proxy=rep(0, 364),
                  weather=weather,
                  ii_init=0, dobserved_true=diff(out$dcc) + 1, ddeaths_true=diff(out$ddeaths) + 1)

stan.code <- drop.stan.model.prior(stan.model.master)

fit <- stan(model_code=stan.code, data=stan.data, control=list(max_treedepth=15))


fit <- tryCatch({
    stan(model_code=stan.code, data=stan.data, open_progress=F, control=list(max_treedepth=15))
}, error=function(e) {
    NULL
})
