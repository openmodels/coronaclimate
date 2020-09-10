setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

format.regtbl <- function(mu, se) {
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

weather <- c('absh', 'r', 't2m', 'tp')
e.ols.mu <- c(0.005168, -0.004100, -0.010323, -0.000976)
e.ols.se <- c(0.003857, 0.000994, 0.005625, 0.000454)

tbl <- data.frame(OLS=c(format.regtbl(e.ols.mu[1], e.ols.se[1]),
                        format.regtbl(e.ols.mu[2], e.ols.se[2]),
                        format.regtbl(e.ols.mu[3], e.ols.se[3]),
                        format.regtbl(e.ols.mu[4], e.ols.se[4]), rep(NA, 8)))

for (model in c('0817-run1', '0817', '0817-noprior', '0907')) {
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

names(tbl) <- c("OLS", "With Prior 1", "With Prior 2", "Without Prior", "Omega Effect")
row.names(tbl) <- c(paste0('e.', c(weather[1], " ", weather[2], "  ", weather[3], "   ", weather[4], "    ")),
                    paste0('o.', c(weather[1], " ", weather[2], "  ", weather[3], "   ", weather[4], "    ")))
xtable(tbl)
