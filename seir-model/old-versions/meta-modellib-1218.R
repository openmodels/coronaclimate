wt.mean <- function(x,wt) {
	s = which(is.finite(x*wt)); wt = wt[s]; x = x[s] #remove NA info
	return( sum(wt * x)/sum(wt) ) #return the mean
}

wt.var <- function(x,wt) {
	s = which(is.finite(x + wt)); wt = wt[s]; x = x[s] #remove NA info
	xbar = wt.mean(x,wt) #get the weighted mean
	return( sum(wt *(x-xbar)^2)*(sum(wt)/(sum(wt)^2-sum(wt^2))) ) #return the variance
}

wt.sd <- function(x,wt) {
	return( sqrt(wt.var(x,wt)) ) #return the standard deviation
}

## Geospatial meta-analysis

bounds <- list("alpha"=c(0, 10), "invgamma"=c(2, 100), "invsigma"=c(2, 100),
               "invkappa"=c(1, 100), "invtheta"=c(1, 100),
       	       "omega"=c(0, 1), "mobility_slope"=c(-1, 10),
	       "deathrate"=c(0, 1), "deathomegaplus"=c(0, 1), "deathlearning"=c(-.01, 0), "error"=c(0, 10),
               "logbeta"=c(-10, 0), "logomega"=c(-10, 0), "eein"=c(0, 1e4),
               "e.absh"=c(-.1, .1), "e.r"=c(-.1, .1), "e.t2m"=c(-.1, .1), "e.tp"=c(-.1, .1), "e.ssrd"=c(-.1, .1), "e.utci"=c(-.1, .1),
	       "o.absh"=c(-.1, .1), "o.r"=c(-.1, .1), "o.t2m"=c(-.1, .1), "o.tp"=c(-.1, .1), "o.ssrd"=c(-.1, .1), "o.utci"=c(-.1, .1))

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

    sigma <- subdfx$sd
    sigma[!is.na(subdfx$rhat)] <- subdfx$sd * subdfx$rhat[!is.na(subdfx$rhat)]
    if (weight == 'pop') {
        stan.data <- list(I=nrow(subdfx), beta=subdfx$mu, sigma=subdfx$sd, weight=subdfx$population / mean(subdfx$population),
                          lobound=min(subdfx$mu), hibound=max(subdfx$mu),
                          mu_prior=weighted.mean(subdfx$mu, subdfx$population / subdfx$sd^2),
                          mu_prior_sd=sqrt(wt.var(subdfx$mu, subdfx$population) + weighted.mean(subdfx$sd, subdfx$population)^2) / sqrt(nrow(subdfx)),
                          tau_prior=wt.sd(subdfx$mu, subdfx$population))
    } else if (weight == 'nobs') {
        stan.data <- list(I=nrow(subdfx), beta=subdfx$mu, sigma=subdfx$sd, weight=subdfx$nobs / mean(subdfx$nobs),
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
        recorded.base <- subdfx[, c('regid', 'param', 'mu', 'sd', 'ci2.5', 'ci25', 'ci50', 'ci75', 'ci97.5', 'Country', 'Region', 'Locality', 'population', 'nobs', 'lowest_level', 'implausible', 'rhat')]
        recorded.base$group <- "Raw"
	recorded.glob <- data.frame(regid=paste(country, region, ""), param, mu=NA, sd=NA, ci2.5=NA, ci25=NA, ci50=NA, ci75=NA, ci97.5=NA, Country=country, Region=region, Locality="", population=sum(subdfx$population), nobs=sum(subdfx$nobs), lowest_level=0, implausible=max(subdfx$implausible), rhat=NA, group="Combined")
        return(rbind(recorded.base, recorded.glob))
    }

    thetas <- la0$theta
    thetas[thetas < bounds[[param]][1] & thetas > bounds[[param]][2]] <- NA
    for (cc in 1:ncol(thetas)) {
        if (sum(!is.na(thetas[, cc])) == 0)
            thetas[, cc] <- seq(bounds[[param]][1], bounds[[param]][2], length.out=nrow(thetas))
    }

    rhats <- stan_rhat(fit0)$data

    subdfx$metamu <- apply(thetas, 2, mean)
    subdfx$metasd <- apply(thetas, 2, sd)
    subdfx$metaci2.5 <- apply(thetas, 2, function(x) quantile(x, .025))
    subdfx$metaci25 <- apply(thetas, 2, function(x) quantile(x, .25))
    subdfx$metaci50 <- apply(thetas, 2, function(x) quantile(x, .5))
    subdfx$metaci75 <- apply(thetas, 2, function(x) quantile(x, .75))
    subdfx$metaci97.5 <- apply(thetas, 2, function(x) quantile(x, .975))

    ## Add on global
    recorded.base <- subdfx[, c('regid', 'param', 'mu', 'sd', 'ci2.5', 'ci25', 'ci50', 'ci75', 'ci97.5', 'Country', 'Region', 'Locality', 'population', 'nobs', 'lowest_level', 'implausible', 'rhat')]
    recorded.base$group <- "Raw"
    recorded.meta <- subdfx[, c('regid', 'param', 'metamu', 'metasd', 'metaci2.5', 'metaci25', 'metaci50', 'metaci75', 'metaci97.5', 'Country', 'Region', 'Locality', 'population', 'nobs', 'lowest_level', 'implausible')]
    names(recorded.meta) <- c('regid', 'param', 'mu', 'sd', 'ci2.5', 'ci25', 'ci50', 'ci75', 'ci97.5', 'Country', 'Region', 'Locality', 'population', 'nobs', 'lowest_level', 'implausible')
    recorded.meta$rhat <- mean(rhats[, 1])
    recorded.meta$group <- "Combined"
    recorded.glob <- data.frame(regid=paste(country, region, ""), param, mu=mean(la0$mu), sd=sd(la0$mu), ci2.5=quantile(la0$mu, .025), ci25=quantile(la0$mu, .25), ci50=quantile(la0$mu, .5), ci75=quantile(la0$mu, .75), ci97.5=quantile(la0$mu, .975), Country=country, Region=region, Locality="", population=sum(subdfx$population), nobs=sum(subdfx$nobs), lowest_level=0, implausible=max(subdfx$implausible), rhat=NA, group="Combined")
    recorded <- rbind(recorded.base, recorded.meta, recorded.glob)

    recorded
}

## allrecorded <- read.csv(outfile)
allrecorded <- data.frame()
for (param in unique(results$param[!is.na(results$mu)])) {
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
            subdf4 <- subset(subdf3, Region == region & Locality != "")
            recorded.region <- estimate.region(subdf4, param, country, region)
            allrecorded <- rbind(allrecorded, recorded.region[-nrow(recorded.region),])
            write.csv(allrecorded, outfile, row.names=F)
            meta.sublocalities <- rbind(meta.sublocalities, recorded.region[nrow(recorded.region),])
        }

        subdfx <- rbind(subset(subdf3, !(Region %in% has.sublocalities) & Region != ""),
                        meta.sublocalities[, -ncol(meta.sublocalities)]) # drop group
        if (nrow(subdfx) > 1) {
            recorded.country <- estimate.region(subdfx, param, country, "")
            allrecorded <- rbind(allrecorded, recorded.country[-nrow(recorded.country),])
            meta.subregions <- rbind(meta.subregions, recorded.country[nrow(recorded.country),])
        } else {
            allrecorded <- rbind(allrecorded, meta.sublocalities)
        }

        write.csv(allrecorded, outfile, row.names=F)
    }

    recorded.global <- estimate.region(rbind(subset(subdf2, !(Country %in% has.subregions)),
                                              meta.subregions[, -ncol(meta.subregions)]),
                                       param, "", "")
    allrecorded <- rbind(allrecorded, recorded.global)
    write.csv(allrecorded, outfile, row.names=F)
}
