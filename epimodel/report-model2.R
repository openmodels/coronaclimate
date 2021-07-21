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

weather <- c('t2m', 'tp', 'ssrd', 'utci')
e.ols.mu <- c(0.02773106208, -0.03134987114, -0.01027632136, -0.02558036184)
e.ols.se <- c(0.01222135339, 0.008463339282, 0.01140204532, 0.009860867169)

tbl <- data.frame(OLS=c(format.regtbl(e.ols.mu[1], e.ols.se[1]),
                        format.regtbl(e.ols.mu[2], e.ols.se[2]),
                        format.regtbl(e.ols.mu[3], e.ols.se[3]),
                        format.regtbl(e.ols.mu[4], e.ols.se[4]), rep(NA, 8)))

for (model in c('0817-run1', '0817', '0817-noprior', '0907', '0921-pop')) {
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

tbl <- data.frame(OLS=c(format.regtbl(e.ols.mu[1], e.ols.se[1]),
                        format.regtbl(e.ols.mu[2], e.ols.se[2]),
                        format.regtbl(e.ols.mu[3], e.ols.se[3]),
                        format.regtbl(e.ols.mu[4], e.ols.se[4]), rep(NA, 8)))

for (model in c('0907', '0907-pop', '0907-nobs')) {
    df <- read.csv(paste0('../../results/epimodel-meta-', model, '.csv'))
    subdf <- subset(df, Country == "" & Region == "")

    e.epi.mu <- c(subdf$mu[subdf$param == 'e.absh'], subdf$mu[subdf$param == 'e.r'],  subdf$mu[subdf$param == 'e.t2m'], subdf$mu[subdf$param == 'e.tp'])
    e.epi.sd <- c(subdf$sd[subdf$param == 'e.absh'], subdf$sd[subdf$param == 'e.r'], subdf$sd[subdf$param == 'e.t2m'], subdf$sd[subdf$param == 'e.tp'])

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
}

library(xtable)

names(tbl) <- c("OLS", "By Region", "By Population", "By Observations")
row.names(tbl) <- c(paste0('e.', c(weather[1], " ", weather[2], "  ", weather[3], "   ", weather[4], "    ")),
                    paste0('o.', c(weather[1], " ", weather[2], "  ", weather[3], "   ", weather[4], "    ")))
xtable(tbl)
