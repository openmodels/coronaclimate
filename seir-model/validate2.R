## setwd("~/research/coronavirus/code/seir-model")

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

for (rate in seq(-100, 100, by=10)) {
    print(rate)
    weather <- -cos(2*pi*(0:364)/365) + 0.5*rnorm(365)
    trues <- rate * ((0:364) / 365) + .1 * weather + .5 * rnorm(365)

    stan.code <- "
data {
  int<lower=0> T; // time periods
  vector[T] trues;
  vector[T] weather;
  real alpha_prior;
}
parameters {
  real<lower=0> alpha;
  real beta;
  real<lower=0> error;
  vector[T-1] dpreds;
}
transformed parameters {
  vector[T] preds;
  preds[1] = beta * weather[1];
  preds[2:T] = cumulative_sum(dpreds) + beta * weather[2:T];
}
model {
  dpreds ~ normal(0, alpha);
  preds ~ normal(trues, error);
  alpha ~ normal(alpha_prior, alpha_prior);
}"

    stan.data <- list(T=length(trues), trues=trues, weather=weather, alpha_prior=.395 / 100)

    fit <- stan(model_code=stan.code, data=stan.data)
    la <- extract(fit, permute=T)

    save(la, file=paste0("../../results/validate2-", rate, "d.RData"))
}

dynamics <- list()
for (rate in seq(-100, 100, by=10)) {
    print(rate)
    load(paste0("../../results/validate2-", rate, ".RData"))
    la1 <- la
    load(paste0("../../results/validate2-", rate, "b.RData"))
    la2 <- la
    load(paste0("../../results/validate2-", rate, "c.RData"))
    la3 <- la
    load(paste0("../../results/validate2-", rate, "d.RData"))
    dynamics[[as.character(rate)]] <- list('a'=la1[c('preds', 'beta')], 'b'=la2[c('preds', 'beta')], 'c'=la3[c('preds', 'beta')], 'd'=la[c('preds', 'beta')])
}

df <- data.frame()
for (srate in names(dynamics)) {
    print(srate)
    rate <- as.numeric(srate)
    expected <- rate / 365
    for (batch in c('a', 'b', 'c', 'd')) {
        for (ii in 1:100) {
            mod <- lm(yy ~ day, data=data.frame(day=1:365, yy=dynamics[[srate]][[batch]][['preds']][ii,]))
            df <- rbind(df, data.frame(rate=rate, param=c('slope', 'beta'), value=c(mod$coeff[2], dynamics[[srate]][[batch]][['beta']][ii]), expected=c(expected, .1), slopeexpected=expected))
        }
    }
}

library(dplyr)
df2 <- df %>% group_by(rate, param, slopeexpected, expected) %>% summarize(value.loo=quantile(value, .025), value.hii=quantile(value, .975),
                                                                           value.lo=quantile(value, .25), value.hi=quantile(value, .75))
df2$label <- "Error in beta estimate"
df2$label[df2$param == 'slope'] <- "Error in slope estimate"

ggplot(df2, aes(x=slopeexpected)) +
    facet_wrap(~ label, scales="free") +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=-.395 / 100) + geom_vline(xintercept=.395 / 100) +
    geom_ribbon(aes(ymin=value.lo - expected, ymax=value.hi - expected), alpha=.5) +
    geom_ribbon(aes(ymin=value.loo - expected, ymax=value.hii - expected), alpha=.5) +
    scale_x_continuous("True slope value", expand=c(0, 0)) + theme_bw() +
    ylab("Error in estimate")
ggsave("../../figures/validate2.pdf", width=6.5, height=3)

ggplot(subset(df2, param == 'slope'), aes(x=slopeexpected)) +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=-.395 / 100) + geom_vline(xintercept=.395 / 100) +
    geom_ribbon(aes(ymin=value.lo - expected, ymax=value.hi - expected), alpha=.5) +
    geom_ribbon(aes(ymin=value.loo - expected, ymax=value.hii - expected), alpha=.5) +
    scale_x_continuous("True slope value", expand=c(0, 0)) + theme_bw() +
    ylab("Error in estimate of slope")
ggsave("../../figures/validate2-slope.pdf", width=5, height=3)
