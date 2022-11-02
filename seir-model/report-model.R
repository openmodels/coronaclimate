setwd("~/research/coronavirus/code/seir-model")

mainmodel <- '0314'

ols.prior.length <- 15

source("../configs.R")

library(xtable)

format.regtbl <- function(mu, se) {
    if (is.na(mu))
        return(c(NA, NA))
    mustr <- formatC(mu, digits=5, format="f")
    sestr <- paste0("(", formatC(se, digits=5, format="f"), ")")
    if (abs(mu / se) > 2.575829303549)
        mustr <- paste0(mustr, "***")
    else if (abs(mu / se) > 1.959963984540)
        mustr <- paste0(mustr, "**")
    else if (abs(mu / se) > 1.644853626951)
        mustr <- paste0(mustr, "*")
    c(mustr, sestr)
}

alltbls <- list()
allgdfs <- list()

## for (filepath in paste0("../../results/epimodel-meta-", mainmodel, "-", c('full', "noprior"), "-all-nobs-nodel.csv")) {
for (filepath in paste0("../../results/epimodel-meta-", mainmodel, "-noprior-", c('all', "mobile"), "-nobs-nodel.csv")) {
    df <- read.csv(filepath)
    subdf <- subset(df, Country == "" & Region == "")

    e.epi.mu <- sapply(weather, function(ww) subdf$mu[subdf$param == paste0('e.', ww)])
    e.epi.sd <- sapply(weather, function(ww) subdf$sd[subdf$param == paste0('e.', ww)])
    o.epi.mu <- sapply(weather, function(ww) subdf$mu[subdf$param == paste0('o.', ww)])
    o.epi.sd <- sapply(weather, function(ww) subdf$sd[subdf$param == paste0('o.', ww)])

    allgdfs[[filepath]] <- data.frame(weather=rep(weather, 2), mu=c(e.epi.mu, o.epi.mu),
                                      sd=c(e.epi.sd, o.epi.sd), channel=rep(c('Transmission', 'Detection'), each=length(weather)))

    subdf2 <- subset(df, Country != "" & Region == "" & group == "Raw")

    e.epi.mu2 <- sapply(weather, function(ww) mean(subdf2$mu[subdf2$param == paste0('e.', ww)], na.rm=T))
    e.epi.sd2 <- sapply(weather, function(ww) mean(subdf2$sd[subdf2$param == paste0('e.', ww)], na.rm=T))
    o.epi.mu2 <- sapply(weather, function(ww) mean(subdf2$mu[subdf2$param == paste0('o.', ww)], na.rm=T))
    o.epi.sd2 <- sapply(weather, function(ww) mean(subdf2$sd[subdf2$param == paste0('o.', ww)], na.rm=T))

    subdf3 <- subset(df, Country != "" & Region == "" & group == "Combined")

    e.epi.mu3 <- sapply(weather, function(ww) mean(subdf3$mu[subdf3$param == paste0('e.', ww)], na.rm=T))
    e.epi.sd3 <- sapply(weather, function(ww) mean(subdf3$sd[subdf3$param == paste0('e.', ww)], na.rm=T))
    o.epi.mu3 <- sapply(weather, function(ww) mean(subdf3$mu[subdf3$param == paste0('o.', ww)], na.rm=T))
    o.epi.sd3 <- sapply(weather, function(ww) mean(subdf3$sd[subdf3$param == paste0('o.', ww)], na.rm=T))

    tbl <- data.frame()
    for (kk in 1:length(weather))
        tbl <- rbind(tbl, data.frame(OLS=format.regtbl(ols.priors.mu[kk], ols.priors.se[kk]),
                                     EPI.raw=format.regtbl(e.epi.mu2[kk], e.epi.sd2[kk]),
                                     EPI.country=format.regtbl(e.epi.mu3[kk], e.epi.sd3[kk]),
                                     EPI.hyper=format.regtbl(e.epi.mu[kk], e.epi.sd[kk])))
    for (kk in 1:length(weather))
        tbl <- rbind(tbl, data.frame(OLS=c(NA, NA),
                                     EPI.raw=format.regtbl(o.epi.mu2[kk], o.epi.sd2[kk]),
                                     EPI.country=format.regtbl(o.epi.mu3[kk], o.epi.sd3[kk]),
                                     EPI.hyper=format.regtbl(o.epi.mu[kk], o.epi.sd[kk])))

    names(tbl) <- c("OLS", "Epi (Raw)", "Epi (Country)", "Epi (Hyper)")
    row.names(tbl) <- paste0(rep(c('e.', 'o.'), each=length(weather)*2),
                             rep(rep(weather, each=2), 2),
                             rep(c('', ' SD'), length(weather)*2))
    xtable(tbl)

    alltbls[[filepath]] <- tbl
}

## Produce table
tbl.row.names <- sapply(0:(4*length(weather)-1), function(ii) {
    ifelse(ii < 2*length(weather), ifelse(ii %% 2 == 0, paste("Transmission", weather[1 + ii / 2]), ""),
    ifelse(ii %% 2 == 0, paste("Detection", weather[1 + (ii - 2*length(weather)) / 2]), ""))
    })
## tbl <- cbind(tbl.row.names, alltbls[[paste0("../../results/epimodel-meta-", mainmodel, "-noprior-all-nobs-nodel.csv")]][, -2],
##              alltbls[[paste0("../../results/epimodel-meta-", mainmodel, "-full-all-nobs-nodel.csv")]][, -1:-2])
tbl <- cbind(tbl.row.names, alltbls[[paste0("../../results/epimodel-meta-", mainmodel, "-noprior-all-nobs-nodel.csv")]][, -2],
             alltbls[[paste0("../../results/epimodel-meta-", mainmodel, "-noprior-mobile-nobs-nodel.csv")]][, -1:-2])
names(tbl) <- c("", "OLS", "Bayes (Country)", "Bayes (Hyper)", "Bayes (Country)", "Bayes (Hyper)")

print(xtable(tbl), include.rownames=F)

## Produce bars

## gdf <- rbind(data.frame(weather=rep(weather, 2), mu=c(ols.priors.mu, ols.priors.mu), sd=c(ols.priors.se, ols.priors.se), channel="OLS",
##                         panel=rep(c("No Prior", "OLS Prior"), each=length(weather))),
##              cbind(panel="No Prior", allgdfs[[paste0('../../results/epimodel-meta-', mainmodel, '-noprior-all-nobs-nodel.csv')]]),
##              cbind(panel="OLS Prior", allgdfs[[paste0('../../results/epimodel-meta-', mainmodel, '-full-all-nobs-nodel.csv')]]))

gdf <- rbind(data.frame(weather=rep(weather, 2), mu=c(ols.priors.mu, ols.priors.mu), sd=c(ols.priors.se, ols.priors.se),
                        channel="OLS", panel=rep(c("All Regions", "Mobile-only"), each=length(weather))),
             cbind(panel="All Regions", allgdfs[[paste0('../../results/epimodel-meta-', mainmodel, '-noprior-all-nobs-nodel.csv')]]),
             cbind(panel="Mobile-only", allgdfs[[paste0('../../results/epimodel-meta-', mainmodel, '-noprior-mobile-nobs-nodel.csv')]]))

library(ggplot2)

## gdf$panel <- factor(gdf$panel, levels=c('No Prior', 'OLS Prior'))
gdf$panel <- factor(gdf$panel, levels=c('All Regions', 'Mobile-only'))
gdf$channel <- factor(gdf$channel, levels=c('OLS', 'Detection', 'Transmission'))
gdf$weather <- factor(gdf$weather, levels=c('t2m', 'utci', 'ssrd', 'tp'))

## ggplot(subset(gdf, channel != 'OLS'), aes(weather, mu, fill=panel)) +
##     facet_grid(channel ~ ., scales="free_y") +
##     geom_bar(stat='identity', position='dodge') + geom_errorbar(aes(ymin=mu - 1.96*sd, ymax=mu + 1.96*sd), position=position_dodge(width=.875), width=.5) + theme_bw() + xlab("Weather variable") + ylab("Normalized weather response") +
##     scale_fill_discrete(name=NULL) + theme(legend.position="bottom")
## ggsave(paste0("~/Dropbox/Coronavirus and Climate/figures/ols-compare-", mainmodel, ".pdf"), width=3, height=4.5)

## Show with OLS

gp <- ggplot(gdf, aes(weather, mu, colour=channel)) +
    facet_grid(panel ~ ., scales="free_y") +
    geom_point(stat='identity', position=position_dodge(width=.5)) +
    geom_linerange(aes(ymin=mu - 1.96*sd, ymax=mu + 1.96*sd), position=position_dodge(width=.5)) +
    theme_bw() + xlab("Weather variable") + ylab("Normalized weather response") +
    scale_colour_discrete(name=NULL) + theme(legend.position="bottom")
ggsave(paste0("~/Dropbox/Coronavirus and Climate/figures/ols-compare-", mainmodel, ".pdf"), gp, width=3.2, height=4.5)

## Get effect of 1 unit change

setwd("~/research/coronavirus/code/seir-model")

mainmodel <- '0314'
weather <- c('t2m', 'tp', 'ssrd', 'utci')
casespath <- "../../cases/panel_all.csv"

df <- read.csv(casespath)

weatherscales <- apply(df[, weather], 2, sd)

results <- read.csv(paste0("../../results/epimodel-meta-", mainmodel, "-noprior-all-nobs-nodel.csv"))
subres <- subset(results, Country == "" & Region == "")

e.epi.mu <- sapply(weather, function(ww) subres$mu[subres$param == paste0('e.', ww)] / weatherscales[weather == ww])
e.epi.sd <- sapply(weather, function(ww) subres$sd[subres$param == paste0('e.', ww)] / weatherscales[weather == ww])
o.epi.mu <- sapply(weather, function(ww) subres$mu[subres$param == paste0('o.', ww)] / weatherscales[weather == ww])
o.epi.sd <- sapply(weather, function(ww) subres$sd[subres$param == paste0('o.', ww)] / weatherscales[weather == ww])

data.frame(weather=rep(weather, 2), mu=c(e.epi.mu, o.epi.mu), sd=c(e.epi.sd, o.epi.sd), percent=100 * c(e.epi.mu, o.epi.mu), pm=100 * 1.96 * c(e.epi.sd, o.epi.sd), channel=rep(c('Transmission', 'Detection'), each=length(weather)))

subres3 <- subset(results, Country != "" & Region == "" & group == "Combined")

e.epi.t2m.perc <- 100 * subres3$mu[subres3$param == "e.t2m"] / weatherscales[weather == 't2m']
quantile(e.epi.t2m.perc)

o.epi.t2m.perc <- 100 * subres3$mu[subres3$param == "o.t2m"] / weatherscales[weather == 't2m']
quantile(o.epi.t2m.perc)
