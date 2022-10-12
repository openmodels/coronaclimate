## setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

source("../configs.R")

casespath <- "../../cases/panel-prepped_MLI.RData"
outpath <- "../../results/epimodel-0907.csv"
weather <- c('absh', 'r', 't2m', 'tp')
regfilter <- function(rows) T
do.multiproc <- T

load(casespath)

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

stan.model.nodice <- "
data {
  int<lower=0> T; // time periods
  int<lower=0> N; // population
  int<lower=0> K; // # weather preds

  int<lower=1> early0;
  int<lower=1> early1;
  int<lower=1> late0;
  int<lower=1> late1;

  real alpha_prior;
  real invsigma_prior; // incubation period (1 / rate of becoming symptomatic)
  real invgamma_prior; // infectious period
  real beta0_prior;
  real eein_prior;

  matrix[T, K] weather;
  vector[T-1] dmobility_proxy;

  real ii_init; // usually 0, so provide as known
  real dobserved_true[T-1];
}
parameters {
  // parameters
  real<lower=0> alpha; // variation in beta
  real<lower=0> beta0;
  real<lower=2> invsigma; // below 2 and daily step doesn't work
  real<lower=2> invgamma; // below 2 and daily step doesn't work
  real<lower=0, upper=1> omega; // observation rate

  // effect of weather
  vector<lower=-1, upper=1>[K] effect;
  vector<lower=-1, upper=1>[K] omegaeffect;

  // affine transform for mobility
  real<lower=0, upper=10> mobility_slope;

  // early and late cases
  real<lower=0, upper=1> portion_early;

  // latent variables
  vector<lower=0>[T-1] eein;
  vector[T-1] dlogbeta;

  real<lower=0> error;
}
transformed parameters {
  // latent variables
  vector[T] logbeta;
  vector<lower=0>[T] ss; // susceptible
  vector<lower=0>[T-1] new_ee1; // newly exposed
  vector<lower=0>[T] ee1; // exposed
  vector<lower=0>[T] ee2;
  vector<lower=0>[T] ii1; // infected
  vector<lower=0>[T] ii2;

  vector<lower=0>[T-1] dcc; // confirmed cases

  logbeta[1] = log(beta0);
  logbeta[2:T] = log(beta0) + cumulative_sum(dlogbeta);

  ss[1] = N;
  ee1[1] = ii_init;
  ee2[1] = ii_init;
  ii1[1] = ii_init;
  ii2[1] = ii_init;
  for (tt in 2:T) {
    new_ee1[tt-1] = exp(logbeta[tt-1] + dot_product(weather[tt-1], effect))*ss[tt-1]*(ii1[tt-1] + ii2[tt-1]) / N;
    ss[tt] = ss[tt-1] - new_ee1[tt-1];
    ee1[tt] = ee1[tt-1] + new_ee1[tt-1] - 2*ee1[tt-1]/invsigma + eein[tt-1];
    ee2[tt] = ee2[tt-1] + 2*ee1[tt-1]/invsigma - 2*ee2[tt-1]/invsigma;
    ii1[tt] = ii1[tt-1] + 2*ee2[tt-1]/invsigma - 2*ii1[tt-1]/invgamma;
    ii2[tt] = ii2[tt-1] + 2*ii1[tt-1]/invgamma - 2*ii2[tt-1]/invgamma;

    if (tt - early0 > 1) {
      if (tt - late0 > 1)
        dcc[tt-1] = (omega + dot_product(weather[tt-early0-1], omegaeffect)) * (portion_early * sum(new_ee1[max(1, tt-early1-1):(tt-early0-1)]) / (early1-early0+1) + (1 - portion_early) * sum(new_ee1[max(1, tt-late1-1):(tt-late0-1)]) / (late1-late0+1));
      else
        dcc[tt-1] = (omega + dot_product(weather[tt-early0-1], omegaeffect)) * (portion_early * sum(new_ee1[max(1, tt-early1-1):(tt-early0-1)]) / (early1-early0+1));
    } else
      dcc[tt-1] = 0;
  }
}
model {
  // priors
  alpha ~ normal(alpha_prior, alpha_prior);
  invsigma ~ gamma(2, 2 / invsigma_prior);
  invgamma ~ gamma(2, 2 / invgamma_prior);
  beta0 ~ normal(beta0_prior, 1);

  // hyperparameters
  eein ~ exponential(1 / eein_prior);

  // mobility proxy
  dlogbeta ~ normal(mobility_slope * dmobility_proxy, alpha);

  // model fit; add 1 to match data
  dobserved_true ~ lognormal(log(dcc + 1), error);
}"

stan.model.deaths <- "
data {
  int<lower=0> T; // time periods
  int<lower=0> N; // population
  int<lower=0> K; // # weather preds

  int<lower=1> early0;
  int<lower=1> early1;
  int<lower=1> late0;
  int<lower=1> late1;

  real alpha_prior;
  real invsigma_prior; // incubation period (1 / rate of becoming symptomatic)
  real invgamma_prior; // infectious period
  real beta0_prior;
  real eein_prior;

  matrix[T, K] weather;
  vector[T-1] dmobility_proxy;

  real ii_init; // usually 0, so provide as known
  real dobserved_true[T-1];
  real ddeaths_true[T-1];
}
parameters {
  // parameters
  real<lower=0> alpha; // variation in beta
  real<lower=0> beta0;
  real<lower=2> invsigma; // below 2 and daily step doesn't work
  real<lower=2> invgamma; // below 2 and daily step doesn't work
  real<lower=0, upper=1> omega; // observation rate

  real<lower=0, upper=.1> deathrate; // rate of death
  real<lower=0, upper=1> deathomegaplus; // additional rate of reported deaths

  // effect of weather
  vector<lower=-1, upper=1>[K] effect;
  vector<lower=-1, upper=1>[K] omegaeffect;

  // affine transform for mobility
  real<lower=0, upper=10> mobility_slope;

  // early and late cases
  real<lower=0, upper=1> portion_early;

  // latent variables
  vector<lower=0>[T-1] eein;
  vector[T-1] dlogbeta;

  real<lower=0> error;
}
transformed parameters {
  // latent variables
  vector[T] logbeta;
  vector<lower=0>[T] ss; // susceptible
  vector<lower=0>[T-1] new_ee1; // newly exposed
  vector<lower=0>[T] ee1; // exposed
  vector<lower=0>[T] ee2;
  vector<lower=0>[T] ii1; // infected
  vector<lower=0>[T] ii2;

  vector<lower=0>[T-1] dcc; // confirmed cases
  vector<lower=0>[T-1] ddeaths; // deaths

  logbeta[1] = log(beta0);
  logbeta[2:T] = log(beta0) + cumulative_sum(dlogbeta);

  ss[1] = N;
  ee1[1] = ii_init;
  ee2[1] = ii_init;
  ii1[1] = ii_init;
  ii2[1] = ii_init;
  for (tt in 2:T) {
    new_ee1[tt-1] = exp(logbeta[tt-1] + dot_product(weather[tt-1], effect))*ss[tt-1]*(ii1[tt-1] + ii2[tt-1]) / N;
    ss[tt] = ss[tt-1] - new_ee1[tt-1];
    ee1[tt] = ee1[tt-1] + new_ee1[tt-1] - 2*ee1[tt-1]/invsigma + eein[tt-1];
    ee2[tt] = ee2[tt-1] + 2*ee1[tt-1]/invsigma - 2*ee2[tt-1]/invsigma;
    ii1[tt] = ii1[tt-1] + 2*ee2[tt-1]/invsigma - 2*ii1[tt-1]/invgamma;
    ii2[tt] = ii2[tt-1] + 2*ii1[tt-1]/invgamma - 2*ii2[tt-1]/invgamma;

    if (tt - early0 > 1) {
      if (tt - late0 > 1)
        dcc[tt-1] = (omega + dot_product(weather[tt-early0-1], omegaeffect)) * (portion_early * sum(new_ee1[max(1, tt-early1-1):(tt-early0-1)]) / (early1-early0+1) + (1 - portion_early) * sum(new_ee1[max(1, tt-late1-1):(tt-late0-1)]) / (late1-late0+1));
      else
        dcc[tt-1] = (omega + dot_product(weather[tt-early0-1], omegaeffect)) * (portion_early * sum(new_ee1[max(1, tt-early1-1):(tt-early0-1)]) / (early1-early0+1));
    } else
      dcc[tt-1] = 0;

    ddeaths[tt-1] = (2*ii2[tt-1]/invgamma) * deathrate * (omega + (1 - omega) * deathomegaplus);
  }
}
model {
  // priors
  alpha ~ normal(alpha_prior, alpha_prior);
  invsigma ~ gamma(2, 2 / invsigma_prior);
  invgamma ~ gamma(2, 2 / invgamma_prior);
  beta0 ~ normal(beta0_prior, 1);

  // hyperparameters
  eein ~ exponential(1 / eein_prior);

  // mobility proxy
  dlogbeta ~ normal(mobility_slope * dmobility_proxy, alpha);

  // model fit; add 1 to match data
  dobserved_true ~ lognormal(log(dcc + 1), error);
  ddeaths_true ~ lognormal(log(ddeaths + 1), error);
}"

stan.compiled.deaths <- stan_model(model_code=stan.model.deaths)
stan.compiled.nodice <- stan_model(model_code=stan.model.nodice)

## Check that all correlations are negative and sds are on similar order of magnitude
## mobdf <- data.frame()
## for (regid in unique(df$regid)) {
##     subdf <- df[df$regid == regid,]
##     if (sum(!is.na(subdf$mobility_pca1)) > 3)
##         mobdf <- rbind(mobdf, data.frame(cor=cor(subdf$mobility_pca1, subdf$residential_percent_change_from_baseline, use='complete'), sd=sd(subdf$mobility_pca1, na.rm=T)))
## }

if (!do.multiproc) {
    if (!file.exists(outpath)) {
        results <- data.frame()
    } else {
        results <- read.csv(outpath)
        results <- subset(results, !is.na(mu))
    }
}

for (regid in unique(df$regid)) {
    if (!do.multiproc) {
        if (regid %in% results$regid)
            next
    }

    subdf <- df[df$regid == regid,]
    if (!regfilter(subdf))
        next

    if (do.multiproc) {
        ## Check if region is claimed
        regfile <- paste0(outpath, "-", regid)
        if (file.exists(regfile))
            next
        ## Claim this region
        fileConn <- file(regfile)
        writeLines(as.character(Sys.getpid()), fileConn)
        close(fileConn)
    }

    print(regid)

    subdf$Confirmed[is.na(subdf$Confirmed)] <- 0
    while (sum(subdf$Confirmed[-1] < subdf$Confirmed[-nrow(subdf)]) > 0) {
        bads <- c(F, subdf$Confirmed[-1] < subdf$Confirmed[-nrow(subdf)])
        subdf$Confirmed[bads] <- c(NA, subdf$Confirmed[-nrow(subdf)])[bads]
    }

    dmobility <- diff(subdf$mobility_pca1)
    dmobility[is.na(dmobility)] <- 0 # Given the affine intercept

    stan.data <- list(T=nrow(subdf), N=round(subdf$population[1]), K=length(weather),
                      early0=weather.delayA[1], early1=weather.delayA[2],
                      late0=weather.delayB[1], late1=weather.delayB[2],
                      alpha_prior=.395 / 100, eein_prior=1,
                      invsigma_prior=5.2, invgamma_prior=2.9,
                      beta0_prior=2.5 * 2.9, dmobility_proxy=dmobility,
                      weather=demeanlist(subdf[, weather], list(factor(rep('all', nrow(subdf))))),
                      ii_init=0, dobserved_true=diff(subdf$Confirmed) + 1)

    if (sum(!is.na(subdf$Deaths) & !is.na(subdf$Confirmed)) > 10) {
        subdf$Deaths[is.na(subdf$Deaths)] <- 0
        while (sum(subdf$Deaths[-1] < subdf$Deaths[-nrow(subdf)]) > 0) {
            bads <- c(F, subdf$Deaths[-1] < subdf$Deaths[-nrow(subdf)])
            subdf$Deaths[bads] <- c(NA, subdf$Deaths[-nrow(subdf)])[bads]
        }

        stan.data$ddeaths_true <- diff(subdf$Deaths) + 1

        fit <- tryCatch({
	    sampling(stan.compiled.deaths, data=stan.data, open_progress=F, control=list(max_treedepth=15))
	}, error=function(e) {
	    NULL
	})
    } else {
        fit <- tryCatch({
	    sampling(stan.compiled.nodice, data=stan.data, open_progress=F)
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
                    get.paramdf(regid, 'omega', la$omega, rhats),
                    get.paramdf(regid, 'mobility_slope', la$mobility_slope, rhats),
                    get.paramdf(regid, 'portion_early', la$portion_early, rhats),
                    get.paramdf(regid, 'deathrate', la$deathrate, rhats),
                    get.paramdf(regid, 'deathomegaplus', la$deathomegaplus, rhats),
                    get.paramdf(regid, 'error', la$error, rhats),
                    get.paramdf(regid, 'logbeta', la$logbeta, rhats),
                    get.paramdf(regid, 'eein', la$eein, rhats))
    for (kk in 1:length(weather))
        resrow <- rbind(resrow, get.paramdf(regid, paste0('e.', weather[kk]), la$effect[,kk], rhats, rhatparam=paste0('effect\\[', kk, '\\]')))
    for (kk in 1:length(weather))
        resrow <- rbind(resrow, get.paramdf(regid, paste0('o.', weather[kk]), la$omegaeffect[,kk], rhats, rhatparam=paste0('effect\\[', kk, '\\]')))

    if (!do.multiproc) {
        results <- rbind(results, resrow)
        write.csv(results, outpath, row.names=F)
    } else {
        ## Release claim
        write.csv(resrow, regfile, row.names=F)
    }
}
