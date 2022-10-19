library(stringr)

stan.model.master <- "
data {
  int<lower=0> T; // time periods
  int<lower=0> N; // population
  int<lower=0> K; // # weather preds

  real alpha_prior;
  real invsigma_prior; // incubation period (1 / rate of becoming symptomatic)
  real invgamma_prior; // infectious period
  real invkappa_prior; // exposed to testing
  real invtheta_prior; // testing to reported
  real beta0_prior;
  real eein_prior;

  vector[K] total_prior;
  vector<lower=0>[K] total_prior_sd;

  matrix[T, K] weather;
  vector[T-1] dmobility_proxy;

  real ii_init; // usually 0, so provide as known
  vector[T-1] dobserved_true; # 1 more than diff(confirmed)
  vector[T-1] ddeaths_true; # 1 more than diff(deaths)
}
transformed data {
  real observed_true = sum(dobserved_true - 1); // +1 when input
}
parameters {
  // parameters
  real<lower=0> alpha; // variation in beta
  real<lower=0> beta0;
  real<lower=2, upper=17> invsigma; // below 2 and daily step doesn't work
  real<lower=2, upper=17> invgamma; // below 2 and daily step doesn't work
  real<lower=1, upper=17> invkappa; // below 1 and daily step doesn't work
  real<lower=1, upper=17> invtheta; // below 1 and daily step doesn't work
  real<lower=0, upper=1> omega0; // initial observation rate

  real<lower=0, upper=.1> deathrate; // rate of death
  real<lower=-.01, upper=0> deathlearning; // decreases deathrate
  real<lower=0, upper=1> deathomegaplus; // additional rate of reported deaths

  // effect of weather
  vector<lower=-.1, upper=.1>[K] effect;
  vector<lower=-.1, upper=.1>[K] omegaeffect;
  vector<lower=-.1, upper=.1>[6] doweffect6;
  vector<lower=-.1, upper=.1>[6] dowomegaeffect6;

  // affine transform for mobility
  real<lower=0, upper=10> mobility_slope;

  // latent variables
  vector<lower=0>[T-1] eein;
  vector[T-1] dlogbeta;
  vector<lower=-.1, upper=.1>[T-2] dlogomega;

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
  vector<lower=0>[T] qq; // waiting to be tested
  vector<lower=0>[T] rr; // waiting to be reported

  // transformed variables
  vector[7] doweffect;
  vector[7] dowomegaeffect;

  vector[T-1] logomega;
  vector<upper=1>[T-1] omega;
  vector<lower=0>[T-1] dcc; // confirmed cases
  vector<lower=0>[T-1] ddeaths; // deaths

  doweffect[7] = 0;
  dowomegaeffect[7] = 0;
  for (dd in 1:6) {
    doweffect[dd] = doweffect6[dd];
    dowomegaeffect[dd] = dowomegaeffect6[dd];
  }

  logbeta[1] = log(beta0);
  logbeta[2:T] = log(beta0) + cumulative_sum(dlogbeta);

  logomega[1] = log(omega0);
  logomega[2:T-1] = log(omega0) + cumulative_sum(dlogomega);

  ss[1] = N;
  ee1[1] = ii_init;
  ee2[1] = ii_init;
  ii1[1] = ii_init;
  ii2[1] = ii_init;
  qq[1] = ii_init;
  rr[1] = ii_init;
  for (tt in 2:T) {
    new_ee1[tt-1] = exp(logbeta[tt-1] + doweffect[1 + (tt % 7)] + dot_product(weather[tt-1], effect))*ss[tt-1]*(ii1[tt-1] + ii2[tt-1]) / N;
    ss[tt] = ss[tt-1] - new_ee1[tt-1];
    ee1[tt] = ee1[tt-1] + new_ee1[tt-1] - 2*ee1[tt-1]/invsigma + eein[tt-1];
    ee2[tt] = ee2[tt-1] + 2*ee1[tt-1]/invsigma - 2*ee2[tt-1]/invsigma;
    ii1[tt] = ii1[tt-1] + 2*ee2[tt-1]/invsigma - 2*ii1[tt-1]/invgamma;
    ii2[tt] = ii2[tt-1] + 2*ii1[tt-1]/invgamma - 2*ii2[tt-1]/invgamma;

    qq[tt] = qq[tt-1] + 2*ee1[tt-1]/invsigma - qq[tt-1]/invkappa;

    omega[tt-1] = exp(logomega[tt-1] + dowomegaeffect[1 + (tt % 7)] + dot_product(weather[tt-1], omegaeffect));
    rr[tt] = rr[tt-1] + omega[tt-1] * qq[tt-1]/invkappa - rr[tt-1]/invtheta;

    dcc[tt-1] = rr[tt-1]/invtheta;
    ddeaths[tt-1] = (2*ii2[tt-1]/invgamma) * deathrate * exp(tt * deathlearning) * (omega[tt-1] + (1 - omega[tt-1]) * deathomegaplus);
  }
}
model {
  // priors
  alpha ~ normal(alpha_prior, alpha_prior);
  invsigma ~ gamma(2, 2 / invsigma_prior);
  invgamma ~ gamma(2, 2 / invgamma_prior);
  invkappa ~ normal(invkappa_prior, invkappa_prior);
  invtheta ~ normal(invtheta_prior, invtheta_prior);
  beta0 ~ normal(beta0_prior, 1);

  total_prior ~ normal(effect + omegaeffect, total_prior_sd);

  // facts of the data
  target += uniform_lpdf(sum(omega .* (dobserved_true - 1) / observed_true) | observed_true / N, 1);

  // hyperparameters
  eein ~ exponential(1 / eein_prior);
  dlogomega ~ normal(0, alpha);

  // mobility proxy
  dlogbeta ~ normal(mobility_slope * dmobility_proxy, alpha);

  // model fit; add 1 to match data
  dobserved_true ~ lognormal(log(dcc + 1), error);
  ddeaths_true ~ lognormal(log(ddeaths + 1), error);
}"

get.stan.model.deaths <- function() stan.model.master
get.stan.model.nodice <- function() {
    one <- str_replace(str_replace(str_replace(str_replace(stan.model.master, fixed("vector[T-1] ddeaths_true;"), ""),
                                               fixed("real<lower=0, upper=.1> deathrate; // rate of death"), ""),
                                   fixed("real<lower=-.01, upper=0> deathlearning; // decreases deathrate"), ""),
                       fixed("real<lower=0, upper=1> deathomegaplus; // additional rate of reported deaths"), "")
    str_replace(str_replace(str_replace(one, fixed("vector<lower=0>[T-1] ddeaths; // deaths"), ""),
                            fixed("ddeaths[tt-1] = (2*ii2[tt-1]/invgamma) * deathrate * exp(tt * deathlearning) * (omega[tt-1] + (1 - omega[tt-1]) * deathomegaplus);"), ""),
                fixed("ddeaths_true ~ lognormal(log(ddeaths + 1), error);"), "")
}

drop.stan.model.prior <- function(model) {
    str_replace(str_replace(str_replace(model, fixed("vector[K] total_prior;"), ""),
                            fixed("vector<lower=0>[K] total_prior_sd;"), ""),
                fixed("total_prior ~ normal(effect + omegaeffect, total_prior_sd);"), "")
}

drop.stan.model.weather <- function(model) {
    model2 <- drop.stan.model.prior(model)
    one <- str_replace(str_replace(str_replace(str_replace(str_replace(model2, fixed("int<lower=0> K; // # weather preds"), ""),
                                                           fixed("matrix[T, K] weather;"), ""),
                                               fixed("// effect of weather"), ""),
                                   fixed("vector<lower=-.1, upper=.1>[K] effect;"), ""),
                       fixed("vector<lower=-.1, upper=.1>[K] omegaeffect;"), "")
    str_replace(str_replace(one, fixed(" + dot_product(weather[tt-1], effect)"), ""),
                fixed(" + dot_product(weather[tt-1], omegaeffect)"), "")
}

drop.stan.model.omega <- function(model) {
    str_replace(str_replace(str_replace(model, fixed("vector<lower=-.1, upper=.1>[K] omegaeffect;"), ""),
                            fixed(" + dot_product(weather[tt-1], omegaeffect)"), ""),
                fixed(" + omegaeffect"), "")
}
                
drop.stan.model.dlogomega <- function(model) {
    str_replace(str_replace(str_replace(model, fixed("vector<lower=-.1, upper=.1>[T-2] dlogomega;"), ""),
                            fixed("logomega[2:T-1] = log(omega0) + cumulative_sum(dlogomega);"), "for (tt in 2:T-1) logomega[tt] = log(omega0);"),
                fixed("dlogomega ~ normal(0, alpha);"), "")
}
