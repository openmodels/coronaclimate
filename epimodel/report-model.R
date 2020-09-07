setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

weather <- c('absh', 'r', 't2m', 'tp')
e.ols.mu <- c(0.005168, -0.004100, -0.010323, -0.000976)
e.ols.se <- c(0.003857, 0.000994, 0.005625, 0.000454)

df <- read.csv('../../results/epimodel-meta-0817.csv')
subdf <- subset(df, Country == "" & Region == "")

e.epi.mu <- c(subdf$mu[subdf$param == 'e.absh'], subdf$mu[subdf$param == 'e.r'],  subdf$mu[subdf$param == 'e.t2m'], subdf$mu[subdf$param == 'e.tp'])
e.epi.sd <- c(subdf$sd[subdf$param == 'e.absh'], subdf$sd[subdf$param == 'e.r'], subdf$sd[subdf$param == 'e.t2m'], subdf$sd[subdf$param == 'e.tp'])

subdf2 <- subset(df, Country != "" & Region == "" & group == "Raw")
e.epi.mu2 <- c(mean(subdf2$mu[subdf2$param == 'e.absh']), mean(subdf2$mu[subdf2$param == 'e.r']), mean(subdf2$mu[subdf2$param == 'e.t2m']), mean(subdf2$mu[subdf2$param == 'e.tp']))
e.epi.sd2 <- c(mean(subdf2$sd[subdf2$param == 'e.absh']), mean(subdf2$sd[subdf2$param == 'e.r']), mean(subdf2$sd[subdf2$param == 'e.t2m']), mean(subdf2$sd[subdf2$param == 'e.tp']))

subdf3 <- subset(df, Country != "" & Region == "" & group == "Combined")
e.epi.mu3 <- c(mean(subdf3$mu[subdf3$param == 'e.absh']), mean(subdf3$mu[subdf3$param == 'e.r']), mean(subdf3$mu[subdf3$param == 'e.t2m']), mean(subdf3$mu[subdf3$param == 'e.tp']))
e.epi.sd3 <- c(mean(subdf3$sd[subdf3$param == 'e.absh']), mean(subdf3$sd[subdf3$param == 'e.r']), mean(subdf3$sd[subdf3$param == 'e.t2m']), mean(subdf3$sd[subdf3$param == 'e.tp']))


library(xtable)

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

tbl <- data.frame(OLS=c(format.regtbl(e.ols.mu[1], e.ols.se[1]),
                        format.regtbl(e.ols.mu[2], e.ols.se[2]),
                        format.regtbl(e.ols.mu[3], e.ols.se[3]),
                        format.regtbl(e.ols.mu[4], e.ols.se[4])),
                  EPI.raw=c(format.regtbl(e.epi.mu2[1], e.epi.sd2[1]),
                         format.regtbl(e.epi.mu2[2], e.epi.sd2[2]),
                         format.regtbl(e.epi.mu2[3], e.epi.sd2[3]),
                         format.regtbl(e.epi.mu2[4], e.epi.sd2[4])),
                  EPI.country=c(format.regtbl(e.epi.mu3[1], e.epi.sd3[1]),
                         format.regtbl(e.epi.mu3[2], e.epi.sd3[2]),
                         format.regtbl(e.epi.mu3[3], e.epi.sd3[3]),
                         format.regtbl(e.epi.mu3[4], e.epi.sd3[4])),
                  EPI.hyper=c(format.regtbl(e.epi.mu[1], e.epi.sd[1]),
                        format.regtbl(e.epi.mu[2], e.epi.sd[2]),
                        format.regtbl(e.epi.mu[3], e.epi.sd[3]),
                        format.regtbl(e.epi.mu[4], e.epi.sd[4])))

names(tbl) <- c("OLS", "Epi (Raw)", "Epi (Country)", "Epi (Hyper)")
row.names(tbl) <- c(weather[1], " ", weather[2], "  ", weather[3], "   ", weather[4], "    ")
xtable(tbl)
