## setwd("~/research/coronavirus/code/seir-model")

source("../configs.R")

source("forward-0314.R")

version <- "0314"

source(paste0("modellib-", version, ".R"))

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

params <- list(invsigma=5.2, invgamma=2.9, invkappa=7, invtheta=7,
               beta0=2.5 * 2.9, doweffect6=rep(0, 6), dowomegaeffect6=rep(0, 6),
               logbeta=rate * ((0:364) / 365) - 0.7, logomega=(-rate / 5) * ((0:364) / 365) - 0.7,
               effect=0.02, omegaeffect=-0.02, deathrate=.05, deathlearning=-.005, deathomegaplus=.2,
               eein=rep(0, 365))

## Generate simulated data

for (rate in c(-.25, .25)) { #seq(-1, 1, by=.5)

weather <- matrix(-cos(2*pi*(0:364)/365) + 0.5*rnorm(365), 365, 1)

data <- list(T=365, N=1e6, K=1, weather=weather, ii_init=1)

out <- forward(data, params)

plot(out$TT, out$dcc)

## Generate simulated data

stan.data <- list(T=365, N=1e6, K=1,
                  alpha_prior=.395 / 100, eein_prior=(22 / 30) / 11e6 * 1e6,
                  invsigma_prior=5.2, invgamma_prior=2.9,
                  invkappa_prior=7, invtheta_prior=7,
                  beta0_prior=2.5 * 2.9, dmobility_proxy=rep(0, 364),
                  weather=weather,
                  ii_init=0, dobserved_true=out$dcc[-1] + 1, ddeaths_true=out$ddeaths[-1] + 1)

stan.code <- drop.stan.model.prior(stan.model.master)

fit <- stan(model_code=stan.code, data=stan.data, control=list(max_treedepth=15))
la <- extract(fit, permute=T)

save(la, file=paste0("../../results/validate-", rate, ".RData"))
}

allla <- list()
dynamics <- list()
for (rate in c(-1, -.5, -.25, 0, .25, .5, 1)) {
    print(rate)
    load(paste0("../../results/validate-", rate, ".RData"))
    if (is.null(allla)) {
        for (name in names(la)) {
            if (length(dim(la[[name]])) == 1 || name %in% c('effect', 'omegaeffect'))
                allla[[name]] <- as.numeric(la[[name]])
        }
        dynamics[[as.character(rate)]] <- la[c('logbeta', 'logomega', 'dcc', 'ddeaths')]
    } else {
        for (name in names(la)) {
            if (length(dim(la[[name]])) == 1 || name %in% c('effect', 'omegaeffect'))
                allla[[name]] <- c(allla[[name]], la[[name]])
        }
        dynamics[[as.character(rate)]] <- la[c('logbeta', 'logomega', 'dcc', 'ddeaths')]
    }
}

pdfs <- data.frame()
trues <- data.frame()
for (param in c('invsigma', 'invgamma', 'invkappa', 'invtheta',
                'effect', 'omegaeffect', 'deathrate', 'deathlearning',
                'deathomegaplus')) {
    pdfs <- rbind(pdfs, data.frame(param, value=as.numeric(la[[param]])))
    trues <- rbind(trues, data.frame(param, value=as.numeric(params[[param]])))
}

ggplot(pdfs, aes(value)) +
    facet_wrap(~ param, scales='free') + geom_histogram() +
    geom_vline(data=trues, aes(xintercept=value)) + theme_bw()

## 1 la

ribbons <- data.frame()
lines <- data.frame()
for (param in c('logbeta', 'logomega', 'dcc', 'ddeaths')) {
    ribbons <- rbind(ribbons, data.frame(param, TT=1:dim(la[[param]])[2], cilo=apply(la[[param]], 2, function(xx) quantile(xx, .25)), cihi=apply(la[[param]], 2, function(xx) quantile(xx, .75)), mu=apply(la[[param]], 2, mean)))
    if (param %in% c('logbeta', 'logomega'))
        lines <- rbind(lines, data.frame(param, TT=1:length(params[[param]]), yy=params[[param]]))
    else
        lines <- rbind(lines, data.frame(param, TT=1:length(out[[param]]), yy=out[[param]]))
}

ggplot(ribbons, aes(TT)) +
    facet_wrap(~ param, scales='free') +
    geom_ribbon(aes(ymin=cilo, ymax=cihi), alpha=.5) + geom_line(aes(y=mu, colour='predicted')) +
    geom_line(data=lines, aes(y=yy, colour='true')) +
    theme_bw()

## Use full dynamics

## Compare slopes

df <- data.frame()
for (srate in names(dynamics)) {
    print(srate)
    rate <- as.numeric(srate)
    for (param in c('logbeta', 'logomega')) {
        expected <- ifelse(param == 'logbeta', rate / 365, (-rate / 5) / 365)
        nday <- ifelse(param == 'logbeta', 365, 364)
        for (ii in 1:400) { # 4000
            mod <- lm(yy ~ day, data=data.frame(day=1:nday, yy=dynamics[[srate]][[param]][ii,]))
            df <- rbind(df, data.frame(rate=rate, param, slope=mod$coeff[2], expected))
        }
    }
}

ggplot(df, aes(x=slope)) +
    facet_grid(~ param + rate, scales='free') +
    geom_histogram() + geom_vline(aes(xintercept=expected))

library(dplyr)
df2 <- df %>% group_by(param, rate, expected) %>% summarize(slope.lo=quantile(slope, .25), slope.hi=quantile(slope, .75))

ggplot(df2, aes(x=rate)) +
    facet_wrap(~ param) +
    geom_line(aes(y=expected)) +
    geom_ribbon(aes(ymin=slope.lo, ymax=slope.hi), alpha=.5)
