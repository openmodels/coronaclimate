setwd("~/Dropbox/Coronavirus and Climate")

library(dplyr)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

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

load("cases/panel-prepped.RData")

get.paramdf <- function(regid, param, lax) {
    data.frame(regid, param, mu=mean(lax), sd=sd(lax), ci2.5=quantile(lax, .025), ci25=quantile(lax, .25),
               ci50=quantile(lax, .50), ci75=quantile(lax, .75), ci97.5=quantile(lax, .975))
}

results <- data.frame()
for (regid in unique(df$regid)) {
    print(regid)

    subdf <- df[df$regid == regid,]
    subdf$Confirmed[is.na(subdf$Confirmed)] <- 0

    stan.data <- list(T=nrow(subdf), N=subdf$population[1], K=3,
                      alpha_prior=.395, eein_prior=10,
                      invsigma_prior=5.2, invgamma_prior=2.9, invkappa_prior=6.1,
                      beta0_prior=2.5 * 2.9,
                      weather=demeanlist(subdf[, c('absh', 'r', 'tp')], list(factor(rep('all', nrow(subdf))))),
                      ii_init=0, dobserved_true=diff(subdf$Confirmed))

    fit <- stan(model_code=stan.model4, data=stan.data, open_progress=F)

    la <- extract(fit, permute=T)

    resrow <- rbind(get.paramdf(regid, 'alpha', la$alpha),
                    get.paramdf(regid, 'invsigma', la$invsigma),
                    get.paramdf(regid, 'invgamma', la$invgamma),
                    get.paramdf(regid, 'invkappa', la$invkappa),
                    get.paramdf(regid, 'omega', la$omega),
                    get.paramdf(regid, 'e.absh', la$effect[,1]),
                    get.paramdf(regid, 'e.r', la$effect[,2]),
                    get.paramdf(regid, 'e.tp', la$effect[,3]),
                    get.paramdf(regid, 'error', la$error))
    results <- rbind(results, resrow)
}
