setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

library(dplyr)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

stan.model1 <- "
data {
  int<lower=0> T;
  int<lower=0> N;

  real alpha_prior;
  real invsigma_prior; // incubation period (1 / rate of becoming symptomatic)
  real invgamma_prior; // infectious period
  real invkappa_prior; // delay from onset to report
  real beta0_prior;
  real fout_prior;
  real eein_prior;

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
  vector<lower=0>[T] ee1;
  vector<lower=0>[T] ee2;
  vector<lower=0>[T] ii1;
  vector<lower=0>[T] ii2;
  vector<lower=0>[T] qq;
  vector<lower=0>[T-1] dcc;

  logbeta[1] = log(beta0);
  ss[1] = N;
  ee1[1] = ii_init;
  ee2[1] = ii_init;
  ii1[1] = ii_init;
  ii2[1] = ii_init;
  qq[1] = ii_init;
  for (tt in 2:T) {
    logbeta[tt] = logbeta[tt-1] + dlogbeta[tt-1];
    ss[tt] = ss[tt-1] - exp(logbeta[tt-1])*ss[tt-1]*(ii1[tt-1] + ii2[tt-1]) / N;
    ee1[tt] = ee1[tt-1] + (1 - fout[tt-1]) * exp(logbeta[tt-1])*ss[tt-1]*(ii1[tt-1] + ii2[tt-1])/N - 2*ee1[tt-1]/invsigma + eein[tt-1];
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

df <- read.csv("hubei_confirmed_cases.csv")
df2 <- df %>% group_by(date) %>% summarize(total.case=sum(total_case))
df2$date <- as.Date(df2$date)

stan.data <- list(T=10 + nrow(df2), N=11e6, alpha_prior=.395),
                  fout_prior=3300/11e6, eein_prior=10, ii_init=1,
                  invsigma_prior=5.2, invgamma_prior=2.9, invkappa_prior=6.1,
                  beta0_prior=2.5 * 2.9, dobserved_true=c(rep(0, 10), diff(df2$total.case)))

fit <- stan(model_code=stan.model1, data=stan.data)

la <- extract(fit, permute=T)

get.byt <- function(name, date, draws) {
    mu <- apply(draws, 2, mean)
    qlo <- apply(draws, 2, function(x) quantile(x, .25))
    qhi <- apply(draws, 2, function(x) quantile(x, .75))
    data.frame(name, date, mu, qlo, qhi)
}

alldates <- c(seq(df2$date[1] - 10, df2$date[1] - 1, by=1), df2$date)

pdf <- rbind(get.byt('exposed', alldates, la$ee1 + la$ee2),
             get.byt('infected', alldates, la$ii1 + la$ii2),
             get.byt('total', alldates[-1], t(apply(la$qq[, -1] - la$qq[, -ncol(la$qq)] + la$dcc, 1, cumsum))),
             get.byt('confirmed', alldates[-1], t(apply(la$dcc * matrix(la$omega, 4000, dim(la$dcc)[2]), 1, cumsum))))

library(ggplot2)
ggplot(pdf, aes(date, mu, colour=name)) +
    geom_line() + geom_ribbon(data=pdf[pdf$name == 'confirmed',], aes(ymin=qlo, ymax=qhi), alpha=.5) +
    geom_point(data=df2, aes(y=total.case, colour='reported')) +
    scale_x_date(expand=c(0, 0)) + theme_bw() + scale_y_log10(limits=c(1, max(pdf$exposed))) +
    xlab(NULL) + ylab(NULL) + scale_colour_discrete(name=NULL)

pdf <- rbind(get.byt('total', alldates[-1], la$qq[, -1] - la$qq[, -ncol(la$qq)] + la$dcc),
             get.byt('confirmed', alldates[-1], la$dcc * matrix(la$omega, 4000, dim(la$dcc)[2])))

ggplot(pdf, aes(date, mu, colour=name)) +
    geom_line() + geom_ribbon(data=pdf[pdf$name == 'confirmed',], aes(ymin=qlo, ymax=qhi), alpha=.5) +
    geom_point(data=data.frame(date=alldates[-1], total.case=c(rep(0, 10), diff(df2$total.case))), aes(y=total.case, colour='reported')) +
    scale_x_date(expand=c(0, 0)) + theme_bw() + scale_y_log10(limits=c(1, max(pdf$total))) +
    xlab(NULL) + ylab(NULL) + scale_colour_discrete(name=NULL)

## Compare priors and results
pdf <- data.frame(x=c(la$beta0 / la$invgamma, rnorm(4000, stan.data$beta0_prior, 1) / rgamma(4000, 2, 2 / stan.data$invgamma_prior),
                      la$alpha, rnorm(4000, stan.data$alpha_prior, .5), la$invsigma, rgamma(4000, 2, 2 / stan.data$invsigma_prior),
                      la$invgamma, rgamma(4000, 2, 2 / stan.data$invgamma_prior),
                      la$invkappa, rexp(4000, 1 / stan.data$invkappa_prior),
                      la$omega, runif(4000, 0, 1)),
                  parameter=rep(c("Initial R0", "R0 Variance", "Incubation Period",
                                  "Infectious Period", "Reporting Delay", "Reporting Fraction"), each=8000),
                  src=rep(c(rep("Fit", 4000), rep("Prior", 4000)), 6))
pdf$x[pdf$parameter == "Initial R0" & pdf$x > 30] <- NA
pdf$x[pdf$parameter == "Reporting Delay" & pdf$x > 40] <- NA
pdf$x[pdf$parameter == "R0 Variance" & pdf$x < 0] <- NA

ggplot(pdf, aes(x, colour=src)) +
    facet_wrap(~ parameter, scales="free") + geom_density() + theme_bw() +
    xlab(NULL)

