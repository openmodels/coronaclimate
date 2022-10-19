casespath <- "../../cases/panel-prepped.RData"
outpath <- "../../results/epimodel-0728.csv"
weather <- c('absh', 'r', 'tp')
regfilter <- function(rows) T
do.multiproc <- T

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

  real ii_init; // usually 0, so provide as known
  real dobserved_true[T-1];
}
parameters {
  // parameters
  real<lower=0> alpha; // variation in beta and omega
  real<lower=0> beta0;
  real<lower=2> invsigma; // below 2 and daily step doesn't work
  real<lower=2> invgamma; // below 2 and daily step doesn't work
  real<lower=0, upper=1> omega0; // initial observation rate

  // effect of weather
  vector<lower=-20, upper=20>[K] effect;

  // early and late cases
  real<lower=0, upper=1> portion_early;

  // latent variables
  vector<lower=0>[T-1] eein;
  vector[T-1] dlogbeta;
  vector[T-2] dlogomega;

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

  vector<lower=0>[T-early0] qqearly; // symptomatic cases able to be reported early
  vector<lower=0>[T-early0-1] dccearly; // cases confirmed early
  vector<lower=0>[T-late0] qqlate; // symptomatic cases able to be reported late
  vector<lower=0>[T-late0-1] dcclate; // cases confirmed late
  vector<lower=0>[T-1] dcc; // confirmed cases

  vector[T-1] logomega;
  vector[T-1] omega;

  logbeta[1] = log(beta0);
  logbeta[2:T] = log(beta0) + cumulative_sum(dlogbeta);

  logomega[1] = log(omega0);
  logomega[2:T-1] = log(omega0) + cumulative_sum(dlogomega);

  omega = exp(logomega) ./ (1 + exp(logomega));

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

    dcc[tt-1] = 0; // initialize to 0
  }

  qqearly[1] = 0;
  for (ss in 2:(T-early0)) {
    // tt = ss + early0 - 1; new_ee1[1] corresponds to tt = 2
    dccearly[ss-1] = qqearly[ss-1] * omega * portion_early;
    qqearly[ss] = qqearly[ss-1] + new_ee1[ss+early0-2] - (qqearly[ss-1] - dccearly[ss-1])/(early1-early0+1);
    dcc[ss+early0-2] += dccearly[ss-1];
  }
  qqlate[1] = 0;
  for (ss in 2:(T-late0)) {
    // ss=2, tt=late0+1, want to get people exposed in tt=1, who then waited in qqearly for early1-early0+1;
    // so in period tt want to grab ssearly = (tt-early0) + early1-early0+1 =
    // tt = sslate + late0 - 1 => grab at ssearly = (sslate+late0-early0-1) + early1-early0+1
    dcclate[ss-1] = qqlate[ss-1] * omega * (1 - portion_early);
    qqlate[ss] = qqlate[ss-1] + (qqearly[(ss+late0-early0-1) + early1-early0+1] - dccearly[(ss+late0-early0-1) + early1-early0+1])/(early1-early0+1) - (qqlate[ss-1] - dcclate[ss-1])/(late1-late0+1);
    dcc[ss+late00-2] += dcclate[ss-1];
  }
}
model {
  // priors
  alpha ~ normal(alpha_prior, .5);
  invsigma ~ gamma(2, 2 / invsigma_prior);
  invgamma ~ gamma(2, 2 / invgamma_prior);
  beta0 ~ normal(beta0_prior, 1);

  // hyperparameters
  eein ~ exponential(1 / eein_prior);
  dlogbeta ~ normal(0, alpha);
  dlogomega ~ normal(0, alpha);

  // model fit
  dobserved_true ~ lognormal(dcc, error);
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
  real invkappa_prior; // delay from onset to death
  real beta0_prior;
  real eein_prior;

  matrix[T, K] weather;

  real ii_init; // usually 0, so provide as known
  real dobserved_true[T-1];
  real ddeaths_true[T-1];
}
parameters {
  // parameters
  real<lower=0> alpha; // variation in beta and omega
  real<lower=0> beta0;
  real<lower=2> invsigma; // below 2 and daily step doesn't work
  real<lower=2> invgamma; // below 2 and daily step doesn't work
  real<lower=1> invkappa; // below 1 and daily step doesn't work
  real<lower=0, upper=1> omega0; // initial observation rate

  real<lower=0, upper=.1> deathrate; // rate of death
  real<lower=0, upper=1> deathomegaplus; // additional rate of reported deaths

  // effect of weather
  vector<lower=-20, upper=20>[K] effect;

  // early and late cases
  real<lower=0, upper=1> portion_early;

  // latent variables
  vector<lower=0>[T-1] eein;
  vector[T-1] dlogbeta;
  vector[T-2] dlogomega;

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

  vector<lower=0>[T-early0] qqearly; // symptomatic cases able to be reported early
  vector<lower=0>[T-early0-1] dccearly; // cases confirmed early
  vector<lower=0>[T-late0] qqlate; // symptomatic cases able to be reported late
  vector<lower=0>[T-late0-1] dcclate; // cases confirmed late
  vector<lower=0>[T-1] dcc; // confirmed cases

  vector[T-1] logomega;
  vector[T-1] omega;
  vector<lower=0>[T-1] dcc; // confirmed cases
  vector<lower=0>[T-1] ddeaths; // deaths

  logbeta[1] = log(beta0);
  logbeta[2:T] = log(beta0) + cumulative_sum(dlogbeta);

  logomega[1] = log(omega0);
  logomega[2:T-1] = log(omega0) + cumulative_sum(dlogomega);

  omega = exp(logomega) ./ (1 + exp(logomega));

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

    dcc[tt-1] = 0; // initialize to 0
  }

  qqearly[1] = 0;
  for (ss in 2:(T-early0)) {
    // tt = ss + early0 - 1; new_ee1[1] corresponds to tt = 2
    dccearly[ss-1] = qqearly[ss-1] * omega * portion_early;
    qqearly[ss] = qqearly[ss-1] + new_ee1[ss+early0-2] - (qqearly[ss-1] - dccearly[ss-1])/(early1-early0+1);
    dcc[ss+early0-2] += dccearly[ss-1];
  }
  qqlate[1] = 0;
  for (ss in 2:(T-late0)) {
    // ss=2, tt=late0+1, want to get people exposed in tt=1, who then waited in qqearly for early1-early0+1;
    // so in period tt want to grab ssearly = (tt-early0) + early1-early0+1 =
    // tt = sslate + late0 - 1 => grab at ssearly = (sslate+late0-early0-1) + early1-early0+1
    dcclate[ss-1] = qqlate[ss-1] * omega * (1 - portion_early);
    qqlate[ss] = qqlate[ss-1] + (qqearly[(ss+late0-early0-1) + early1-early0+1] - dccearly[(ss+late0-early0-1) + early1-early0+1])/(early1-early0+1) - (qqlate[ss-1] - dcclate[ss-1])/(late1-late0+1);
    dcc[ss+late00-2] += dcclate[ss-1];
  }

  ddeaths = dcc .* deathrate;
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
  dlogomega ~ normal(0, alpha);

  // model fit
  dobserved_true ~ lognormal(dcc, error);
  ddeaths_true ~ lognormal(ddeaths .* (omega + (1 - omega) * deathomegaplus), error);
}"


load(casespath)

get.paramdf <- function(regid, param, lax) {
    if (is.null(lax) || length(lax) == 0)
        return(data.frame(regid, param, mu=NA, sd=NA, ci2.5=NA, ci25=NA, ci50=NA, ci75=NA, ci97.5=NA))
    data.frame(regid, param, mu=mean(lax), sd=sd(lax), ci2.5=quantile(lax, .025), ci25=quantile(lax, .25),
               ci50=quantile(lax, .50), ci75=quantile(lax, .75), ci97.5=quantile(lax, .975))
}

stan.compiled.deaths <- stan_model(model_code=stan.model.deaths)
stan.compiled.nodice <- stan_model(model_code=stan.model.nodice)

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

    stan.data <- list(T=nrow(subdf), N=round(subdf$population[1]), K=length(weather),
                      alpha_prior=.395 / 100, eein_prior=1,
                      invsigma_prior=5.2, invgamma_prior=2.9, invkappa_prior=6.1,
                      beta0_prior=2.5 * 2.9,
                      weather=demeanlist(subdf[, weather], list(factor(rep('all', nrow(subdf))))),
                      ii_init=0, dobserved_true=diff(subdf$Confirmed) + 1)

    if (sum(!is.na(subdf$Deaths) & !is.na(subdf$Confirmed)) > 10) {
        subdf$Deaths[is.na(subdf$Deaths)] <- 0
        while (sum(subdf$Deaths[-1] < subdf$Deaths[-nrow(subdf)]) > 0) {
            bads <- c(F, subdf$Deaths[-1] < subdf$Deaths[-nrow(subdf)])
            subdf$Deaths[bads] <- c(NA, subdf$Deaths[-nrow(subdf)])[bads]
        }

        stan.data$ddeaths_true <- diff(subdf$Deaths) + 1

        fit <- sampling(stan.compiled.deaths, data=stan.data, open_progress=F)
    } else {
        fit <- sampling(stan.compiled.nodice, data=stan.data, open_progress=F)
    }

    la <- extract(fit, permute=T)

    resrow <- rbind(get.paramdf(regid, 'alpha', la$alpha),
                    get.paramdf(regid, 'invsigma', la$invsigma),
                    get.paramdf(regid, 'invgamma', la$invgamma),
                    get.paramdf(regid, 'invkappa', la$invkappa),
                    get.paramdf(regid, 'omega', la$omega),
                    get.paramdf(regid, 'deathrate', la$deathrate),
                    get.paramdf(regid, 'deathomegaplus', la$deathomegaplus),
                    get.paramdf(regid, 'e.absh', la$effect[,1]),
                    get.paramdf(regid, 'e.r', la$effect[,2]),
                    get.paramdf(regid, 'e.tp', la$effect[,3]),
                    get.paramdf(regid, 'error', la$error))

    if (!do.multiproc) {
        results <- rbind(results, resrow)
        write.csv(results, outpath, row.names=F)
    } else {
        ## Release claim
        write.csv(resrow, regfile, row.names=F)
    }
}
