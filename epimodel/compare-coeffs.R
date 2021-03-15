setwd("~/research/coronavirus")

df <- read.csv("results-saved/epimodel-meta-0105noprior-all-nobs.csv")
df.glo <- subset(df, Country == "" & Region == "")
df.glo$weather <- NA
df.glo$weather[grep("e\\.", df.glo$param)] <- gsub("e\\.", "", df.glo$param[grep("e\\.", df.glo$param)])
df.glo$weather[grep("o\\.", df.glo$param)] <- gsub("o\\.", "", df.glo$param[grep("o\\.", df.glo$param)])
df.glo$channel <- NA
df.glo$channel[grep("e\\.", df.glo$param)] <- "beta"
df.glo$channel[grep("o\\.", df.glo$param)] <- "omega"

coeffs <- rbind(read.csv("code/epimodel/model-0210newprior-two.csv"),
                read.csv("code/epimodel/model-0210newprior-one.csv"))
coeffs$weather <- sapply(coeffs$X, function(x) strsplit(as.character(x), "_")[[1]][1])
coeffs$channel <- NA
coeffs$channel[grep("_L1", coeffs$X)] <- "omega"
coeffs$channel[grep("_L2", coeffs$X)] <- "both"

library(dplyr)

df.omega <- coeffs %>% inner_join(df.glo)
df.both <- subset(coeffs, channel == 'both') %>% left_join(df.glo, by='weather', suffix=c('.ols', '.epi'))
df.omega$channel.epi <- df.omega$channel
names(df.omega)[names(df.omega) == 'channel'] <- 'channel.ols'

df.compare <- rbind(df.omega, df.both)

library(ggplot2)

ggplot(df.compare, aes(mu, Estimate, colour=channel.epi)) +
    facet_wrap(~ channel.ols) +
    geom_point(aes(shape=weather)) +
    #geom_segment(aes(y=Estimate - Cluster.s.e., yend=Estimate + Cluster.s.e., xend=mu)) +
    #geom_segment(aes(x=mu - sd, xend=mu + sd, yend=Estimate)) +
    #geom_smooth(method='lm', formula=y ~ 0 + x, se=F, fullrange=T) +
    theme_bw() + xlab("Epidemiological estimates (beta and omega)") +
    ylab("OLS estimates (omega and both)")

## Find out what we would expect to see from these params

source("code/epimodel/forward-0105.R")

params <- list()
for (param in unique(df.glo$param))
    params[[param]] <- df.glo$mu[df.glo$param == param]

params[['logbeta']] <- rep(params[['logbeta']], 365)
params[['logomega']] <- rep(params[['logomega']], 365)
params[['doweffect6']] <- rep(0, 6)
params[['dowomegaeffect6']] <- rep(0, 6)
params[['eein']] <- rep(0, 364)

results <- data.frame()
for (weatvar in c('absh', 'ssrd', 't2m', 'tp', 'utci')) {
    param.e <- paste0('e.', weatvar)
    param.o <- paste0('o.', weatvar)
    params[['effect']] <- params[[param.e]]
    params[['omegaeffect']] <- params[[param.o]]

    weather <- matrix(0, 365, 1)
    data <- list(T=365, N=1, K=1, weather=weather, ii_init=1 / 10e6)

    baseline <- get.dlog(data, params)
    
    impresp <- get.impresp(baseline, 100, params)
    results <- rbind(results, data.frame(weather=weatvar, time=1:364, impresp=impresp))
}

res.avg <- results %>% group_by(weather) %>%
    summarize(simulated=c(mean(impresp[time - 100 > 1.5 & time - 100 < 3.5]),
                          mean(impresp[time - 100 > 3.5 & time - 100 < 12.5]),
                          mean(impresp[time - 100 > 12.5 & time - 100 < 30])),
              channel=c('omega', 'both', 'three'))

df.sim <- coeffs %>% left_join(res.avg)
df.sim$channel <- factor(df.sim$channel, levels=c('omega', 'both'))

ggplot(df.sim[!is.na(df.sim$channel),], aes(simulated, Estimate)) +
    facet_wrap(~ channel) +
    geom_point(aes(shape=weather)) +
    #geom_segment(aes(y=Estimate - Cluster.s.e., yend=Estimate + Cluster.s.e., xend=mu)) +
    #geom_segment(aes(x=mu - sd, xend=mu + sd, yend=Estimate)) +
    geom_abline(slope=1) +
    theme_bw() + xlab("Epidemiological estimates under simulation") +
    ylab("OLS estimates")
