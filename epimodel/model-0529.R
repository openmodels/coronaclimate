setwd("~/Dropbox/Coronavirus and Climate")

library(dplyr)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

load("cases/panel-prepped.RData")

poldf <- read.csv("policy/policyml.csv")
df$fex <- poldf$fex

subdf <- subset(df, regid == "United Kingdom  ")
subdf$Confirmed[is.na(subdf$Confirmed)] <- 0

stan.model2 <- "
data {
  int<lower=0> T; // time periods
  int<lower=0> N; // population
  int<lower=0> K; // # weather preds
  int<lower=0> L; // # policy regimes

  real alpha_prior;
  real invsigma_prior; // incubation period (1 / rate of becoming symptomatic)
  real invgamma_prior; // infectious period
  real invkappa_prior; // delay from onset to report
  real beta0_prior;
  real fout_prior;
  real eein_prior;

  matrix[T, K] weather;
  int<lower=1, upper=L> regime[T];

  real ii_init; // usually 0, so provide as known
  real dobserved_true[T-1];
}
parameters {
  // parameters
  real<lower=0> alpha;
  real<lower=0> beta0;
  real<lower=2> invsigma; // below 2 and daily step doesn't work
  real<lower=2> invgamma; // below 2 and daily step doesn't work
  real<lower=1> invkappa; // below 1 and daily step doesn't work
  real<lower=0, upper=1> omega;

  // effect of weather
  vector[K] effect;
  vector[L] fex;

  // latent variables
  vector<lower=0, upper=1>[T-1] fout;
  vector<lower=0>[T-1] eein;
  vector[T-1] dlogbeta;

  real<lower=0> error;
}
transformed parameters {
  // latent variables
  vector[T] logbeta;
  vector<lower=0>[T] ss;
  vector<lower=0>[T-1] new_ee1;
  vector<lower=0>[T] ee1;
  vector<lower=0>[T] ee2;
  vector<lower=0>[T] ii1;
  vector<lower=0>[T] ii2;
  vector<lower=0>[T] qq;
  vector<lower=0>[T-1] dcc;

  logbeta[1] = log(beta0);
  logbeta[2:T] = log(beta0) + cumulative_sum(dlogbeta);

  ss[1] = N;
  ee1[1] = ii_init;
  ee2[1] = ii_init;
  ii1[1] = ii_init;
  ii2[1] = ii_init;
  qq[1] = ii_init;
  for (tt in 2:T) {
    //logbeta[tt] = logbeta[tt-1] + dlogbeta[tt-1];
    new_ee1[tt-1] = exp(logbeta[tt-1] + dot_product(weather[tt-1], effect) + fex[regime[tt-1]])*ss[tt-1]*(ii1[tt-1] + ii2[tt-1]) / N;
    ss[tt] = ss[tt-1] - new_ee1[tt-1];
    ee1[tt] = ee1[tt-1] + (1 - fout[tt-1]) * new_ee1[tt-1] - 2*ee1[tt-1]/invsigma + eein[tt-1];
    ee2[tt] = ee2[tt-1] + 2*ee1[tt-1]/invsigma - 2*ee2[tt-1]/invsigma;
    ii1[tt] = ii1[tt-1] + 2*ee2[tt-1]/invsigma - 2*ii1[tt-1]/invgamma;
    ii2[tt] = ii2[tt-1] + 2*ii1[tt-1]/invgamma - 2*ii2[tt-1]/invgamma;
    qq[tt] = qq[tt-1] + 2*ee2[tt-1]*exp(-1/(invgamma*invkappa))/invsigma - qq[tt-1]/invkappa;
    dcc[tt-1] = qq[tt-1]/invkappa;
  }
}
model {
  // priors
  alpha ~ normal(alpha_prior, .5);
  invsigma ~ gamma(2, 2 / invsigma_prior);
  invgamma ~ gamma(2, 2 / invgamma_prior);
  invkappa ~ exponential(1 / invkappa_prior);
  beta0 ~ normal(beta0_prior, 1);

  // hyperparameters
  fout ~ exponential(1 / fout_prior);
  eein ~ exponential(1 / eein_prior);
  dlogbeta ~ normal(0, alpha);

  // model fit
  dobserved_true ~ normal(dcc * omega, error);
}"

stan.data <- list(T=nrow(subdf), N=subdf$population[1], K=3, L=length(unique(subdf$fex)),
                  alpha_prior=.395, fout_prior=3300/11e6, eein_prior=10,
                  invsigma_prior=5.2, invgamma_prior=2.9, invkappa_prior=6.1,
                  beta0_prior=2.5 * 2.9,
                  weather=subdf[, c('absh', 'r', 'tp')],
                  regime=as.numeric(factor(subdf$fex, levels=unique(subdf$fex))),
                  ii_init=0,
                  dobserved_true=diff(subdf$Confirmed))

fit <- stan(model_code=stan.model2, data=stan.data)

stan.data <- list(T=nrow(subdf), N=subdf$population[1], K=3, L=1,
                  alpha_prior=.395, fout_prior=3300/11e6, eein_prior=10,
                  invsigma_prior=5.2, invgamma_prior=2.9, invkappa_prior=6.1,
                  beta0_prior=2.5 * 2.9,
                  weather=subdf[, c('absh', 'r', 'tp')],
                  regime=rep(1, nrow(subdf)),
                  ii_init=0,
                  dobserved_true=diff(subdf$Confirmed))

fit <- stan(model_code=stan.model2, data=stan.data)

stan.model3 <- "
data {
  int<lower=0> T; // time periods
  int<lower=0> N; // population
  int<lower=0> K; // # weather preds

  real alpha_prior;
  real invsigma_prior; // incubation period (1 / rate of becoming symptomatic)
  real invgamma_prior; // infectious period
  real invkappa_prior; // delay from onset to report
  real beta0_prior;
  real fout_prior;
  real eein_prior;

  matrix[T, K] weather;

  real ii_init; // usually 0, so provide as known
  real dobserved_true[T-1];
}
parameters {
  // parameters
  real<lower=0> alpha;
  real<lower=0> beta0;
  real<lower=2> invsigma; // below 2 and daily step doesn't work
  real<lower=2> invgamma; // below 2 and daily step doesn't work
  real<lower=1> invkappa; // below 1 and daily step doesn't work
  real<lower=0, upper=1> omega;

  // effect of weather
  vector[K] effect;

  // latent variables
  vector<lower=0, upper=1>[T-1] fout;
  vector<lower=0>[T-1] eein;
  vector[T-1] dlogbeta;

  real<lower=0> error;
}
transformed parameters {
  // latent variables
  vector[T] logbeta;
  vector<lower=0>[T] ss;
  vector<lower=0>[T-1] new_ee1;
  vector<lower=0>[T] ee1;
  vector<lower=0>[T] ee2;
  vector<lower=0>[T] ii1;
  vector<lower=0>[T] ii2;
  vector<lower=0>[T] qq;
  vector<lower=0>[T-1] dcc;

  logbeta[1] = log(beta0);
  logbeta[2:T] = log(beta0) + cumulative_sum(dlogbeta);

  ss[1] = N;
  ee1[1] = ii_init;
  ee2[1] = ii_init;
  ii1[1] = ii_init;
  ii2[1] = ii_init;
  qq[1] = ii_init;
  for (tt in 2:T) {
    //logbeta[tt] = logbeta[tt-1] + dlogbeta[tt-1];
    new_ee1[tt-1] = exp(logbeta[tt-1] + dot_product(weather[tt-1], effect))*ss[tt-1]*(ii1[tt-1] + ii2[tt-1]) / N;
    ss[tt] = ss[tt-1] - new_ee1[tt-1];
    ee1[tt] = ee1[tt-1] + (1 - fout[tt-1]) * new_ee1[tt-1] - 2*ee1[tt-1]/invsigma + eein[tt-1];
    ee2[tt] = ee2[tt-1] + 2*ee1[tt-1]/invsigma - 2*ee2[tt-1]/invsigma;
    ii1[tt] = ii1[tt-1] + 2*ee2[tt-1]/invsigma - 2*ii1[tt-1]/invgamma;
    ii2[tt] = ii2[tt-1] + 2*ii1[tt-1]/invgamma - 2*ii2[tt-1]/invgamma;
    qq[tt] = qq[tt-1] + 2*ee2[tt-1]*exp(-1/(invgamma*invkappa))/invsigma - qq[tt-1]/invkappa;
    dcc[tt-1] = qq[tt-1]/invkappa;
  }
}
model {
  // priors
  alpha ~ normal(alpha_prior, .5);
  invsigma ~ gamma(2, 2 / invsigma_prior);
  invgamma ~ gamma(2, 2 / invgamma_prior);
  invkappa ~ exponential(1 / invkappa_prior);
  beta0 ~ normal(beta0_prior, 1);

  // hyperparameters
  fout ~ exponential(1 / fout_prior);
  eein ~ exponential(1 / eein_prior);
  dlogbeta ~ normal(0, alpha);

  // model fit
  dobserved_true ~ normal(dcc * omega, error);
}"

library(lfe)
stan.data <- list(T=nrow(subdf), N=subdf$population[1], K=3,
                  alpha_prior=.395, fout_prior=3300/11e6, eein_prior=10,
                  invsigma_prior=5.2, invgamma_prior=2.9, invkappa_prior=6.1,
                  beta0_prior=2.5 * 2.9,
                  weather=demeanlist(subdf[, c('absh', 'r', 'tp')], list(factor(rep('all', nrow(subdf))))),
                  ii_init=0, dobserved_true=diff(subdf$Confirmed))

fit <- stan(model_code=stan.model3, data=stan.data)

stan.data <- list(T=nrow(subdf), N=subdf$population[1], K=0,
                  alpha_prior=.395, fout_prior=3300/11e6, eein_prior=10,
                  invsigma_prior=5.2, invgamma_prior=2.9, invkappa_prior=6.1,
                  beta0_prior=2.5 * 2.9,
                  weather=demeanlist(subdf[, c()], list(factor(rep('all', nrow(subdf))))),
                  ii_init=0, dobserved_true=diff(subdf$Confirmed))

fit <- stan(model_code=stan.model3, data=stan.data)

stan.model4 <- "
data {
  int<lower=0> T; // time periods
  int<lower=0> N; // population
  int<lower=0> K; // # weather preds

  real alpha_prior;
  real invsigma_prior; // incubation period (1 / rate of becoming symptomatic)
  real invgamma_prior; // infectious period
  real invkappa_prior; // delay from onset to report
  real beta0_prior;
  real eein_prior;

  matrix[T, K] weather;

  real ii_init; // usually 0, so provide as known
  real dobserved_true[T-1];
}
parameters {
  // parameters
  real<lower=0> alpha;
  real<lower=0> beta0;
  real<lower=2> invsigma; // below 2 and daily step doesn't work
  real<lower=2> invgamma; // below 2 and daily step doesn't work
  real<lower=1> invkappa; // below 1 and daily step doesn't work
  real<lower=0, upper=1> omega;

  // effect of weather
  vector[K] effect;

  // latent variables
  vector<lower=0>[T-1] eein;
  vector[T-1] dlogbeta;

  real<lower=0> error;
}
transformed parameters {
  // latent variables
  vector[T] logbeta;
  vector<lower=0>[T] ss;
  vector<lower=0>[T-1] new_ee1;
  vector<lower=0>[T] ee1;
  vector<lower=0>[T] ee2;
  vector<lower=0>[T] ii1;
  vector<lower=0>[T] ii2;
  vector<lower=0>[T] qq;
  vector<lower=0>[T-1] dcc;

  logbeta[1] = log(beta0);
  logbeta[2:T] = log(beta0) + cumulative_sum(dlogbeta);

  ss[1] = N;
  ee1[1] = ii_init;
  ee2[1] = ii_init;
  ii1[1] = ii_init;
  ii2[1] = ii_init;
  qq[1] = ii_init;
  for (tt in 2:T) {
    new_ee1[tt-1] = exp(logbeta[tt-1] + dot_product(weather[tt-1], effect))*ss[tt-1]*(ii1[tt-1] + ii2[tt-1]) / N;
    ss[tt] = ss[tt-1] - new_ee1[tt-1];
    ee1[tt] = ee1[tt-1] + new_ee1[tt-1] - 2*ee1[tt-1]/invsigma + eein[tt-1];
    ee2[tt] = ee2[tt-1] + 2*ee1[tt-1]/invsigma - 2*ee2[tt-1]/invsigma;
    ii1[tt] = ii1[tt-1] + 2*ee2[tt-1]/invsigma - 2*ii1[tt-1]/invgamma;
    ii2[tt] = ii2[tt-1] + 2*ii1[tt-1]/invgamma - 2*ii2[tt-1]/invgamma;
    qq[tt] = qq[tt-1] + 2*ee2[tt-1]*exp(-1/(invgamma*invkappa))/invsigma - qq[tt-1]/invkappa;
    dcc[tt-1] = qq[tt-1]/invkappa;
  }
}
model {
  // priors
  alpha ~ normal(alpha_prior, .5);
  invsigma ~ gamma(2, 2 / invsigma_prior);
  invgamma ~ gamma(2, 2 / invgamma_prior);
  invkappa ~ exponential(1 / invkappa_prior);
  beta0 ~ normal(beta0_prior, 1);

  // hyperparameters
  eein ~ exponential(1 / eein_prior);
  dlogbeta ~ normal(0, alpha);

  // model fit
  dobserved_true ~ normal(dcc * omega, error);
}"

stan.data <- list(T=nrow(subdf), N=subdf$population[1], K=3,
                  alpha_prior=.395, eein_prior=10,
                  invsigma_prior=5.2, invgamma_prior=2.9, invkappa_prior=6.1,
                  beta0_prior=2.5 * 2.9,
                  weather=demeanlist(subdf[, c('absh', 'r', 'tp')], list(factor(rep('all', nrow(subdf))))),
                  ii_init=0, dobserved_true=diff(subdf$Confirmed))

fit <- stan(model_code=stan.model4, data=stan.data)

la <- extract(fit, permute=T)


