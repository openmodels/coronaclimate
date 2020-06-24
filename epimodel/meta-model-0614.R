setwd("~/Dropbox/Coronavirus and Climate")

library(dplyr)
library(ggplot2)
library(PBSmapping)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

results <- read.csv("results/epimodel-0616.csv")

ggplot(results, aes(mu)) +
    facet_wrap(~ param, scales='free') +
    geom_histogram() + xlab(NULL) + ylab(NULL) + theme_bw()

## Geospatial meta-analysis

bounds <- list("omega"=c(0, 1))

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

stan.model <- "
data {
  int<lower=0> I; // number of studies
  real beta[I]; // estimated treatment effects
  real<lower=0> sigma[I]; // s.e. of effect estimates
  real cos_lon[I]; // cos(longitude)
  real sin_lon[I]; // sin(longitude)
  real div_lat[I]; // latitude / 90
  real<lower=0> geo_prior;
  real lobound;
  real hibound;
}
parameters {
  real<lower=lobound, upper=hibound> mu;
  real gamma_cosl;
  real gamma_sinl;
  real gamma_xlat;
  real<lower=0> tau;
  real eta[I];
}
transformed parameters {
  real theta[I];
  for (ii in 1:I)
    theta[ii] = mu * (1 + gamma_cosl * cos_lon[ii] + gamma_sinl * sin_lon[ii] + gamma_xlat * div_lat[ii]) + tau * eta[ii];
}
model {
  target += normal_lpdf(eta | 0, 1);
  target += normal_lpdf(beta | theta, sigma);
  gamma_cosl ~ normal(0, geo_prior);
  gamma_sinl ~ normal(0, geo_prior);
  gamma_xlat ~ normal(0, geo_prior);
}
"

load("cases/panel-prepped.RData")
df2 <- df[, c('regid', 'ALPHA.3', 'population', 'lowest_level', 'implausible')] %>% group_by(regid) %>% summarize(ALPHA.3=ALPHA.3[1], population=mean(population), lowest_level=min(lowest_level), implausible=max(implausible))

shp <- importShapefile("shapefiles/gadm36_levels_simple/adm0.shp")
centroids <- calcCentroid(shp, rollup=1)
polydata <- attr(shp, 'PolyData')
polydata2 <- cbind(polydata, centroids)

df3 <- df2 %>% left_join(polydata2[, -1], by=c('ALPHA.3'='GID_0'))

allrecorded <- data.frame()
for (param in unique(results$param)) {
    subdf <- results[results$param == param & !is.na(results$mu) & !is.na(results$sd),]
    subdf2 <- subdf %>% left_join(df3)

    stan.data <- list(I=nrow(subdf2), beta=subdf2$mu, sigma=subdf2$sd,
                      div_lat=subdf2$Y / 90, cos_lon=cos(subdf2$X * pi / 180),
                      sin_lon=sin(subdf2$Y * pi / 180), geo_prior=0,
                      lobound=bounds[[param]][1], hibound=bounds[[param]][2])

    fit0 <- stan(model_code=stan.model0, data=stan.data,
                 iter=1000, chains=4, control=list(adapt_delta=0.99, max_treedepth=20))
    la0 <- extract(fit0, permute=T)

    stan.data$geo_prior <- sd(la0$mu) # could argue should be sd(la0$theta), but too much freedom

    fit <- stan(model_code=stan.model, data=stan.data,
                iter=1000, chains=4)

    la <- extract(fit, permute=T)

    subdf2$metamu <- apply(la$theta, 2, mean)
    subdf2$metasd <- apply(la$theta, 2, sd)
    subdf2$metaci2.5 <- apply(la$theta, 2, function(x) quantile(x, .025))
    subdf2$metaci25 <- apply(la$theta, 2, function(x) quantile(x, .25))
    subdf2$metaci50 <- apply(la$theta, 2, function(x) quantile(x, .5))
    subdf2$metaci75 <- apply(la$theta, 2, function(x) quantile(x, .75))
    subdf2$metaci97.5 <- apply(la$theta, 2, function(x) quantile(x, .975))

    ## Combine to global level
    subdf3 <- subdf2 %>% group_by(PID) %>% summarize(param=param[1], country=NAME_0[1], mu=mean(mu), sd=mean(sd), ci2.5=mean(ci2.5), ci25=mean(ci25), ci50=mean(ci50), ci75=mean(ci75), ci97.5=mean(ci97.5), metamu=mean(metamu), metasd=mean(metasd), metaci2.5=mean(metaci2.5), metaci25=mean(metaci25), metaci50=mean(metaci50), metaci75=mean(metaci75), metaci97.5=mean(metaci97.5))

    ## Add on global
    recorded.base <- subdf3[, c('param', 'country', 'PID', 'mu', 'sd', 'ci2.5', 'ci25', 'ci50', 'ci75', 'ci97.5')]
    recorded.base$group <- "Raw"
    recorded.meta <- subdf3[, c('param', 'country', 'PID', 'metamu', 'metasd', 'metaci2.5', 'metaci25', 'metaci50', 'metaci75', 'metaci97.5')]
    names(recorded.meta) <- c('param', 'country', 'PID', 'mu', 'sd', 'ci2.5', 'ci25', 'ci50', 'ci75', 'ci97.5')
    recorded.meta$group <- "Combined"
    recorded.glob <- data.frame(param, country="Global", PID=NA, mu=mean(la0$mu), sd=sd(la$mu), ci2.5=quantile(la0$mu, .025), ci25=quantile(la0$mu, .25), ci50=quantile(la0$mu, .5), ci75=quantile(la0$mu, .75), ci97.5=quantile(la0$mu, .975), group="Combined")
    recorded <- rbind(recorded.base, recorded.meta, recorded.glob)

    ggplot(recorded, aes(country, mu, colour=group)) +
        coord_flip() +
        geom_point() + geom_errorbar(aes(ymin=ci25, ymax=ci75)) +
        theme_bw() + theme(axis.text.y=element_text(size=4)) + ylab(param) + xlab(NULL) +
        scale_colour_discrete(name=NULL)

    shp2 <- shp %>% left_join(subdf3[, c('PID', 'metamu', 'metasd')])

    ggplot(shp2, aes(X, Y, fill=metamu, group=paste(PID, SID))) +
        geom_polygon() + scale_y_continuous(name=NULL, limits=c(-60, 85), expand=c(0, 0)) +
        scale_x_continuous(name=NULL, expand=c(0, 0)) + theme_bw() + scale_fill_continuous(name=param)

    allrecorded <- rbind(allrecorded, recorded)
}

ggplot(subset(allrecorded, country == "Global"), aes(param, mu)) +
    coord_flip() +
    geom_point() + geom_errorbar(aes(ymin=ci25, ymax=ci75)) +
    theme_bw() + ylab("Hyper-paramater value and 75% CI") + xlab(NULL)
