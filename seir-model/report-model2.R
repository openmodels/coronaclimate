setwd("~/research/coronavirus/code/epimodel")

format.regtbl <- function(mu, se) {
    if (is.na(mu))
        return(c("NA", "NA"))
    mustr <- formatC(mu, digits=5, format="f")
    if (is.na(se))
        return(c(mustr, "NA"))
    sestr <- paste0("(", formatC(se, digits=5, format="f"), ")")
    if (abs(mu / se) > 2.575829303549)
        mustr <- paste0(mustr, "***")
    else if (abs(mu / se) > 1.959963984540)
        mustr <- paste0(mustr, "**")
    else if (abs(mu / se) > 1.644853626951)
        mustr <- paste0(mustr, "*")
    c(mustr, sestr)
}

source("../configs.R")

tbl <- data.frame(OLS=c(format.regtbl(e.ols.mu[1], e.ols.se[1]),
                        format.regtbl(e.ols.mu[2], e.ols.se[2]),
                        format.regtbl(e.ols.mu[3], e.ols.se[3]),
                        format.regtbl(e.ols.mu[4], e.ols.se[4]), rep(NA, 8)))

for (model in c('0314-noprior', '0314-full3')) {
    df <- read.csv(paste0('../../results/epimodel-meta-', model, '.csv'))
    subdf <- subset(df, Country == "" & Region == "")

    e.epi.mu <- c(subdf$mu[subdf$param == 'e.absh'], subdf$mu[subdf$param == 'e.r'],  subdf$mu[subdf$param == 'e.t2m'], subdf$mu[subdf$param == 'e.tp'])
    e.epi.sd <- c(subdf$sd[subdf$param == 'e.absh'], subdf$sd[subdf$param == 'e.r'], subdf$sd[subdf$param == 'e.t2m'], subdf$sd[subdf$param == 'e.tp'])

    if (model == '0907') {
        o.epi.mu <- c(subdf$mu[subdf$param == 'o.absh'], subdf$mu[subdf$param == 'o.r'],  subdf$mu[subdf$param == 'o.t2m'], subdf$mu[subdf$param == 'o.tp'])
        o.epi.sd <- c(subdf$sd[subdf$param == 'o.absh'], subdf$sd[subdf$param == 'o.r'], subdf$sd[subdf$param == 'o.t2m'], subdf$sd[subdf$param == 'o.tp'])

        tbl[, model] <- c(format.regtbl(e.epi.mu[1], e.epi.sd[1]),
                          format.regtbl(e.epi.mu[2], e.epi.sd[2]),
                          format.regtbl(e.epi.mu[3], e.epi.sd[3]),
                          format.regtbl(e.epi.mu[4], e.epi.sd[4]),
                          format.regtbl(o.epi.mu[1], o.epi.sd[1]),
                          format.regtbl(o.epi.mu[2], o.epi.sd[2]),
                          format.regtbl(o.epi.mu[3], o.epi.sd[3]),
                          format.regtbl(o.epi.mu[4], o.epi.sd[4]))
    } else {
        tbl[, model] <- c(format.regtbl(e.epi.mu[1], e.epi.sd[1]),
                          format.regtbl(e.epi.mu[2], e.epi.sd[2]),
                          format.regtbl(e.epi.mu[3], e.epi.sd[3]),
                          format.regtbl(e.epi.mu[4], e.epi.sd[4]), rep(NA, 8))
    }
}

library(xtable)

names(tbl) <- c("OLS", "With Prior 1", "With Prior 2", "Without Prior", "Omega Effect", "Smooth Omega")
row.names(tbl) <- c(paste0('e.', c(weather[1], " ", weather[2], "  ", weather[3], "   ", weather[4], "    ")),
                    paste0('o.', c(weather[1], " ", weather[2], "  ", weather[3], "   ", weather[4], "    ")))
xtable(tbl)

### Report with different weighting

tbl <- data.frame(OLS=c(format.regtbl(ols.priors.mu[1], ols.priors.se[1]),
                        format.regtbl(ols.priors.mu[2], ols.priors.se[2]),
                        format.regtbl(ols.priors.mu[3], ols.priors.se[3]),
                        format.regtbl(ols.priors.mu[4], ols.priors.se[4]), rep(NA, 8)))

for (model in c('0314-noprior-all-nobs', '0314-noprior-all-pop', '0314-noprior-all-region')) {
    df <- read.csv(paste0('../../results/epimodel-meta-', model, '-nodel.csv'))
    subdf <- subset(df, Country == "" & Region == "")

    e.epi.mu <- c(subdf$mu[subdf$param == 'e.t2m'], subdf$mu[subdf$param == 'e.tp'],  subdf$mu[subdf$param == 'e.ssrd'], subdf$mu[subdf$param == 'e.utci'])
    e.epi.sd <- c(subdf$sd[subdf$param == 'e.t2m'], subdf$sd[subdf$param == 'e.tp'], subdf$sd[subdf$param == 'e.ssrd'], subdf$sd[subdf$param == 'e.utci'])

    o.epi.mu <- c(subdf$mu[subdf$param == 'o.t2m'], subdf$mu[subdf$param == 'o.tp'],  subdf$mu[subdf$param == 'o.ssrd'], subdf$mu[subdf$param == 'o.utci'])
    o.epi.sd <- c(subdf$sd[subdf$param == 'o.t2m'], subdf$sd[subdf$param == 'o.tp'], subdf$sd[subdf$param == 'o.ssrd'], subdf$sd[subdf$param == 'o.utci'])

    tbl[, model] <- c(format.regtbl(e.epi.mu[1], e.epi.sd[1]),
                      format.regtbl(e.epi.mu[2], e.epi.sd[2]),
                      format.regtbl(e.epi.mu[3], e.epi.sd[3]),
                      format.regtbl(e.epi.mu[4], e.epi.sd[4]),
                      format.regtbl(o.epi.mu[1], o.epi.sd[1]),
                      format.regtbl(o.epi.mu[2], o.epi.sd[2]),
                      format.regtbl(o.epi.mu[3], o.epi.sd[3]),
                      format.regtbl(o.epi.mu[4], o.epi.sd[4]))
}

library(xtable)

names(tbl) <- c("OLS", "By Observations", "By Population", "By Region")
row.names(tbl) <- c(paste0('e.', c(weather[1], " ", weather[2], "  ", weather[3], "   ", weather[4], "    ")),
                    paste0('o.', c(weather[1], " ", weather[2], "  ", weather[3], "   ", weather[4], "    ")))
xtable(tbl)
