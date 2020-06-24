## setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

library(dplyr)
library(ggplot2)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

results <- read.csv("../../results/epimodel-0616.csv")
outfile <- "../../results/epimodel-meta-0616.csv"

results$param <- factor(results$param, c('alpha', 'invgamma', 'invkappa', 'invsigma',
                                         'omega', 'deathrate', 'deathomegaplus', 'error',
                                         'e.absh', 'e.r', 'e.tp'))
ggplot(results, aes(mu)) +
    facet_wrap(~ param, scales='free') +
    geom_histogram() + xlab(NULL) + ylab(NULL) + theme_bw()

## Geospatial meta-analysis

bounds <- list("alpha"=c(0, 10), "deathrateplus"=c(0, 1), "deathrate"=c(0, 1),
               "e.absh"=c(-20, 20), "e.r"=c(-20, 20), "e.tp"=c(-20, 20),
               "invgamma"=c(0, 100), "invkappa"=c(0, 100), "invsigma"=c(0, 100), "omega"=c(0, 1))

stan.model0 <- "
data {
  int<lower=0> I; // number of studies
  real beta[I]; // estimated treatment effects
  real<lower=0> sigma[I]; // s.e. of effect estimates
  real lobound;
  real hibound;
}
parameters {
  real<lower=lobound, upper=hibound> mu;
  real<lower=0> tau;
  real eta[I];
}
transformed parameters {
  real theta[I];
  for (ii in 1:I)
    theta[ii] = mu + tau * eta[ii];
}
model {
  target += normal_lpdf(eta | 0, 1);
  target += normal_lpdf(beta | theta, sigma);
}
"

stan.model0.compiled <- stan_model(model_code=stan.model0)

load("../../cases/panel-prepped.RData")
df2 <- df[, c('regid', 'Country', 'Region', 'Locality', 'population', 'lowest_level', 'implausible')] %>% group_by(regid) %>% summarize(Country=Country[1], Region=Region[1], Locality=Locality[1], population=mean(population), lowest_level=min(lowest_level), implausible=max(implausible))

estimate.region <- function(subdfx, param, country, region) {
    print(c(param, country, region))
    stan.data <- list(I=nrow(subdfx), beta=subdfx$mu, sigma=subdfx$sd,
                      lobound=bounds[[param]][1], hibound=bounds[[param]][2])

    fit0 <- sampling(stan.model0.compiled, data=stan.data,
                     iter=1000, chains=4, open_progress=F)
    la0 <- extract(fit0, permute=T)

    thetas <- la0$theta
    thetas[thetas < bounds[[param]][1] & thetas > bounds[[param]][2]] <- NA
    for (cc in 1:ncol(thetas)) {
        if (sum(!is.na(thetas[, cc])) == 0)
            thetas[, cc] <- seq(bounds[[param]][1], bounds[[param]][2], length.out=nrow(thetas))
    }

    subdfx$metamu <- apply(thetas, 2, mean)
    subdfx$metasd <- apply(thetas, 2, sd)
    subdfx$metaci2.5 <- apply(thetas, 2, function(x) quantile(x, .025))
    subdfx$metaci25 <- apply(thetas, 2, function(x) quantile(x, .25))
    subdfx$metaci50 <- apply(thetas, 2, function(x) quantile(x, .5))
    subdfx$metaci75 <- apply(thetas, 2, function(x) quantile(x, .75))
    subdfx$metaci97.5 <- apply(thetas, 2, function(x) quantile(x, .975))

    ## Add on global
    recorded.base <- subdfx[, c('regid', 'param', 'mu', 'sd', 'ci2.5', 'ci25', 'ci50', 'ci75', 'ci97.5', 'Country', 'Region', 'Locality', 'population', 'lowest_level', 'implausible')]
    recorded.base$group <- "Raw"
    recorded.meta <- subdfx[, c('regid', 'param', 'metamu', 'metasd', 'metaci2.5', 'metaci25', 'metaci50', 'metaci75', 'metaci97.5', 'Country', 'Region', 'Locality', 'population', 'lowest_level', 'implausible')]
    names(recorded.meta) <- c('regid', 'param', 'mu', 'sd', 'ci2.5', 'ci25', 'ci50', 'ci75', 'ci97.5', 'Country', 'Region', 'Locality', 'population', 'lowest_level', 'implausible')
    recorded.meta$group <- "Combined"
    recorded.glob <- data.frame(regid=paste(country, region, ""), param, mu=mean(la0$mu), sd=sd(la0$mu), ci2.5=quantile(la0$mu, .025), ci25=quantile(la0$mu, .25), ci50=quantile(la0$mu, .5), ci75=quantile(la0$mu, .75), ci97.5=quantile(la0$mu, .975), Country=country, Region=region, Locality="", population=sum(subdfx$population), lowest_level=0, implausible=max(subdfx$implausible), group="Combined")
    recorded <- rbind(recorded.base, recorded.meta, recorded.glob)

    recorded
}

allrecorded <- data.frame()
for (param in unique(results$param)) {
    if (!(param %in% names(bounds)))
        next

    subdf <- results[results$param == param & !is.na(results$mu) & !is.na(results$sd),]
    subdf2 <- subdf %>% left_join(df2)

    globaldf <- data.frame()
    has.subregions <- names(table(subdf2$Country)[table(subdf2$Country) > 1])
    meta.subregions <- data.frame()
    for (country in has.subregions) {
        subdf3 <- subset(subdf2, Country == country)

        countrydf <- data.frame()
        has.sublocalities <- names(table(subdf3$Region)[table(subdf3$Region) > 1])
        meta.sublocalities <- data.frame()
        for (region in has.sublocalities) {
            subdf4 <- subset(subdf3, Region == region)
            recorded.region <- estimate.region(subdf4, param, country, region)
            allrecorded <- rbind(allrecorded, recorded.region[-nrow(recorded.region),])
            write.csv(allrecorded, outfile, row.names=F)
            meta.sublocalities <- rbind(meta.sublocalities, recorded.region[nrow(recorded.region),])
        }

        subdfx <- rbind(subset(subdf3, !(Region %in% has.sublocalities)),
                        meta.sublocalities[, -ncol(meta.sublocalities)]) # drop group
        if (nrow(subdfx) > 1) {
            recorded.country <- estimate.region(subdfx, param, country, "")
            allrecorded <- rbind(allrecorded, recorded.country[-nrow(recorded.country),])
        } else {
            allrecorded <- rbind(allrecorded, meta.sublocalities)
        }

        write.csv(allrecorded, outfile, row.names=F)

        meta.subregions <- rbind(meta.subregions, recorded.country[nrow(recorded.country),])
    }

    recorded.global <- estimate.region(rbind(subset(subdf2, !(Region %in% has.subregions)),
                                              meta.subregions[, -ncol(meta.subregions)]),
                                       param, "", "")
    allrecorded <- rbind(allrecorded, recorded.global)
    write.csv(allrecorded, outfile, row.names=F)
}
