## setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

library(dplyr)
library(ggplot2)
library(SDMTools)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

weight <- 'pop'

do.display <- F

results <- read.csv("../../results-20200910/epimodel-0907.csv")
outfile <- paste0("../../results-20200910/epimodel-meta-0907-", weight, ".csv")

if (do.display) {
results$param <- factor(results$param, c('alpha', 'invgamma', 'invsigma', 'mobility_slope',
                                         'omega', 'portion_early',
					 'deathrate', 'deathomegaplus', 'error',
					 'logbeta', 'eein',
                                         'e.absh', 'e.r', 'e.t2m', 'e.tp',
					 'o.absh', 'o.r', 'o.t2m', 'o.tp'))

## For plot, drop beyond the 99th
results$showit <- T
for (param in unique(results$param)) {
    limits <- quantile(results$mu[results$param == param], c(.025, .975), na.rm=T)
    results$showit[results$param == param] <- !is.na(results$mu[results$param == param]) & results$mu[results$param == param] > limits[1] & results$mu[results$param == param] < limits[2]
}

results$paramlabel <- as.character(results$param)
results$paramlabel[results$param == 'mobility_slope'] <- "Mobility Adjustment"
results$paramlabel[results$param == 'alpha'] <- "Gradual Adjustment Rate"
results$paramlabel[results$param == 'invsigma'] <- "Incubation Period (days)"
results$paramlabel[results$param == 'invgamma'] <- "Infectious Period (days)"
results$paramlabel[results$param == 'portion_early'] <- "Portion Detected Early"
results$paramlabel[results$param == 'omega'] <- "Recording Rate"
results$paramlabel[results$param == 'deathrate'] <- "Death Rate"
results$paramlabel[results$param == 'deathomegaplus'] <- "Extra Record of Deaths"
results$paramlabel[results$param == 'portion_early'] <- "Portion Reported Early"
results$paramlabel[results$param == 'e.t2m'] <- "Air Temperature Trans."
results$paramlabel[results$param == 'e.tp'] <- "Total Precipitation Trans."
results$paramlabel[results$param == 'e.r'] <- "Relative Humidity Trans."
results$paramlabel[results$param == 'e.absh'] <- "Absolute Humidity Trans."
results$paramlabel[results$param == 'o.t2m'] <- "Air Temperature Detect"
results$paramlabel[results$param == 'o.tp'] <- "Total Precipitation Detect"
results$paramlabel[results$param == 'o.r'] <- "Relative Humidity Detect"
results$paramlabel[results$param == 'o.absh'] <- "Absolute Humidity Detect"

results$paramlabel <- factor(results$paramlabel, levels=c("Gradual Adjustment Rate", "Mobility Adjustment", "Incubation Period (days)", "Infectious Period (days)", "Portion Detected Early", "Recording Rate", "Death Rate", "Extra Record of Deaths", "Portion Reported Early", "Air Temperature Trans.", "Total Precipitation Trans.", "Relative Humidity Trans.", "Absolute Humidity Trans.", "Air Temperature Detect", "Total Precipitation Detect", "Relative Humidity Detect", "Absolute Humidity Detect"))

ggplot(subset(results, param != 'error' & showit), aes(mu)) +
    facet_wrap(~ paramlabel, scales='free') +
    geom_histogram() + xlab(NULL) + ylab(NULL) + theme_bw()

## Scatter plot
library(reshape2)
results2 <- dcast(subset(results, param != 'error' & showit), regid ~ param, value.var='mu')

# Correlation panel
panel.cor <- function(x, y){
    usr <- par("usr"); on.exit(par(usr))
    par(usr=c(0, 1, 0, 1))
    r <- round(cor(x, y, use='complete', method='spearman'), digits=2)
    txt <- as.character(r)
    text(0.5, 0.5, txt)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x, y, pch=19, cex=.01)
}

pairs(results2[,-1], lower.panel=panel.cor, upper.panel=upper.panel)
}

## Geospatial meta-analysis

bounds <- list("alpha"=c(0, 10), "invgamma"=c(2, 100), "invsigma"=c(2, 100),
       	       "omega"=c(0, 1), "mobility_slope"=c(-1, 10), "portion_early"=c(0, 1),
	       "deathrate"=c(0, 1), "deathomegaplus"=c(0, 1), "error"=c(0, 10),
               "e.absh"=c(-1, 1), "e.r"=c(-1, 1), "e.t2m"=c(-1, 1), "e.tp"=c(-1, 1),
	       "o.absh"=c(-1, 1), "o.r"=c(-1, 1), "o.t2m"=c(-1, 1), "o.tp"=c(-1, 1))

stan.model0 <- "
data {
  int<lower=0> I; // number of studies
  vector[I] beta; // estimated treatment effects
  vector<lower=0>[I] sigma; // s.e. of effect estimates
  vector<lower=0>[I] weight; // weight of each study
  real lobound;
  real hibound;
  real mu_prior;
  real<lower=0> mu_prior_sd;
  real<lower=0> tau_prior;
}
parameters {
  real<lower=lobound, upper=hibound> mu;
  real<lower=0> tau;
  vector<lower=lobound, upper=hibound>[I] theta;
}
model {
  for (ii in 1:I)
    target += weight[ii] * normal_lpdf(theta[ii] | mu, tau);
  // theta ~ normal(mu, tau);
  beta ~ normal(theta, sigma);
  mu ~ normal(mu_prior, mu_prior_sd);
  tau ~ cauchy(0, tau_prior);
}
"

stan.model0.compiled <- stan_model(model_code=stan.model0)

stan.data <- list(I=2, beta=c(0, 1), sigma=c(1, 1), weight=c(1, 100),
                  lobound=-100, hibound=100,
                  mu_prior=.5,
                  mu_prior_sd=1000,
                  tau_prior=1000)

fit0 <- sampling(stan.model0.compiled, data=stan.data,
                 iter=2000, chains=4, open_progress=F)
fit0

load("../../cases/panel-prepped_MLI.RData")
df2 <- df[, c('regid', 'Country', 'Region', 'Locality', 'population', 'lowest_level', 'implausible')] %>% group_by(regid) %>% summarize(Country=Country[1], Region=Region[1], Locality=Locality[1], population=mean(population), lowest_level=min(lowest_level), implausible=max(implausible))

estimate.region <- function(subdfx, param, country, region) {
    if (nrow(subdfx) > 100) {
        subsubs <- floor(sqrt(nrow(subdfx)))
        subwhich <- sample(1:subsubs, nrow(subdfx), replace=T)

        recorded <- data.frame()
        subglobs <- data.frame()
        for (ss in 1:subsubs) {
            subrecorded <- estimate.region(subdfx[subwhich == ss,], param, paste0(country, '-', ss), paste(region, '-', ss))
            recorded <- rbind(recorded, subrecorded[-nrow(subrecorded),])
            subglobs <- rbind(subglobs, subrecorded[nrow(subrecorded),])
        }

        fullglob <- estimate.region(subglobs[, -ncol(subglobs)], param, country, region)
        recorded <- rbind(recorded, fullglob)
        return(recorded)
    }

    print(c(param, country, region))
    if (weight == 'pop') {
        stan.data <- list(I=nrow(subdfx), beta=subdfx$mu, sigma=subdfx$sd, weight=subdfx$population,
                          lobound=min(subdfx$mu), hibound=max(subdfx$mu),
                          mu_prior=weighted.mean(subdfx$mu, subdfx$population / subdfx$sd^2),
                          mu_prior_sd=sqrt(wt.var(subdfx$mu, subdfx$population) + weighted.mean(subdfx$sd, subdfx$population)^2) / sqrt(nrow(subdfx)),
                          tau_prior=wt.sd(subdfx$mu, subdfx$population))
    } else {
        stan.data <- list(I=nrow(subdfx), beta=subdfx$mu, sigma=subdfx$sd, weight=rep(1, nrow(subdfx)),
                          lobound=min(subdfx$mu), hibound=max(subdfx$mu),
                          mu_prior=weighted.mean(subdfx$mu, 1 / subdfx$sd^2),
                          mu_prior_sd=sqrt(var(subdfx$mu) + mean(subdfx$sd)^2) / sqrt(nrow(subdfx)),
                          tau_prior=sd(subdfx$mu))
    }

    fit0 <- sampling(stan.model0.compiled, data=stan.data,
                     iter=2000, chains=4, open_progress=F)
    la0 <- extract(fit0, permute=T)
    if (is.null(la0)) {
        recorded.base <- subdfx[, c('regid', 'param', 'mu', 'sd', 'ci2.5', 'ci25', 'ci50', 'ci75', 'ci97.5', 'Country', 'Region', 'Locality', 'population', 'lowest_level', 'implausible', 'rhat')]
        recorded.base$group <- "Raw"
	recorded.glob <- data.frame(regid=paste(country, region, ""), param, mu=NA, sd=NA, ci2.5=NA, ci25=NA, ci50=NA, ci75=NA, ci97.5=NA, Country=country, Region=region, Locality="", population=sum(subdfx$population), lowest_level=0, implausible=max(subdfx$implausible), rhat=NA, group="Combined")
        return(rbind(recorded.base, recorded.glob))
    }

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
    recorded.base <- subdfx[, c('regid', 'param', 'mu', 'sd', 'ci2.5', 'ci25', 'ci50', 'ci75', 'ci97.5', 'Country', 'Region', 'Locality', 'population', 'lowest_level', 'implausible', 'rhat')]
    recorded.base$group <- "Raw"
    recorded.meta <- subdfx[, c('regid', 'param', 'metamu', 'metasd', 'metaci2.5', 'metaci25', 'metaci50', 'metaci75', 'metaci97.5', 'Country', 'Region', 'Locality', 'population', 'lowest_level', 'implausible')]
    names(recorded.meta) <- c('regid', 'param', 'mu', 'sd', 'ci2.5', 'ci25', 'ci50', 'ci75', 'ci97.5', 'Country', 'Region', 'Locality', 'population', 'lowest_level', 'implausible')
    recorded.meta$rhat <- NA
    recorded.meta$group <- "Combined"
    recorded.glob <- data.frame(regid=paste(country, region, ""), param, mu=mean(la0$mu), sd=sd(la0$mu), ci2.5=quantile(la0$mu, .025), ci25=quantile(la0$mu, .25), ci50=quantile(la0$mu, .5), ci75=quantile(la0$mu, .75), ci97.5=quantile(la0$mu, .975), Country=country, Region=region, Locality="", population=sum(subdfx$population), lowest_level=0, implausible=max(subdfx$implausible), rhat=NA, group="Combined")
    recorded <- rbind(recorded.base, recorded.meta, recorded.glob)

    recorded
}

## allrecorded <- read.csv(outfile)
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
