setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

weather <- c('absh', 'ssrd', 't2m', 'tp', 'utci')
e.ols.mu <- c(0.03344, -0.00332, -0.04791, -0.00190, -0.00684)
e.ols.se <- c(0.01493, 0.00423, 0.01812, 0.00135, 0.00632)

df <- read.csv('../../results/epimodel-meta-1030-mobile-pop.csv')
subdf <- subset(df, Country == "" & Region == "")

e.epi.mu <- c(subdf$mu[subdf$param == 'e.absh'], subdf$mu[subdf$param == 'e.ssrd'], subdf$mu[subdf$param == 'e.t2m'], subdf$mu[subdf$param == 'e.tp'], subdf$mu[subdf$param == 'e.utci'])
e.epi.sd <- c(subdf$sd[subdf$param == 'e.absh'], subdf$sd[subdf$param == 'e.ssrd'], subdf$sd[subdf$param == 'e.t2m'], subdf$sd[subdf$param == 'e.tp'], subdf$sd[subdf$param == 'e.utci'])

o.epi.mu <- c(subdf$mu[subdf$param == 'o.absh'], subdf$mu[subdf$param == 'o.ssrd'], subdf$mu[subdf$param == 'o.t2m'], subdf$mu[subdf$param == 'o.tp'], subdf$mu[subdf$param == 'o.utci'])
o.epi.sd <- c(subdf$sd[subdf$param == 'o.absh'], subdf$sd[subdf$param == 'o.ssrd'], subdf$sd[subdf$param == 'o.t2m'], subdf$sd[subdf$param == 'o.tp'], subdf$sd[subdf$param == 'o.utci'])

subdf2 <- subset(df, Country != "" & Region == "" & group == "Raw")
e.epi.mu2 <- c(mean(subdf2$mu[subdf2$param == 'e.absh']), mean(subdf2$mu[subdf2$param == 'e.ssrd']), mean(subdf2$mu[subdf2$param == 'e.t2m']), mean(subdf2$mu[subdf2$param == 'e.tp']), mean(subdf2$mu[subdf2$param == 'e.utci']))
e.epi.sd2 <- c(mean(subdf2$sd[subdf2$param == 'e.absh']), mean(subdf2$sd[subdf2$param == 'e.ssrd']), mean(subdf2$sd[subdf2$param == 'e.t2m']), mean(subdf2$sd[subdf2$param == 'e.tp']), mean(subdf2$sd[subdf2$param == 'e.utci']))

o.epi.mu2 <- c(mean(subdf2$mu[subdf2$param == 'o.absh']), mean(subdf2$mu[subdf2$param == 'o.ssrd']), mean(subdf2$mu[subdf2$param == 'o.t2m']), mean(subdf2$mu[subdf2$param == 'o.tp']), mean(subdf2$mu[subdf2$param == 'o.utci']))
o.epi.sd2 <- c(mean(subdf2$sd[subdf2$param == 'o.absh']), mean(subdf2$sd[subdf2$param == 'o.ssrd']), mean(subdf2$sd[subdf2$param == 'o.t2m']), mean(subdf2$sd[subdf2$param == 'o.tp']), mean(subdf2$sd[subdf2$param == 'o.utci']))

subdf3 <- subset(df, Country != "" & Region == "" & group == "Combined")
e.epi.mu3 <- c(mean(subdf3$mu[subdf3$param == 'e.absh']), mean(subdf3$mu[subdf3$param == 'e.ssrd']), mean(subdf3$mu[subdf3$param == 'e.t2m']), mean(subdf3$mu[subdf3$param == 'e.tp']), mean(subdf3$mu[subdf3$param == 'e.utci']))
e.epi.sd3 <- c(mean(subdf3$sd[subdf3$param == 'e.absh']), mean(subdf3$sd[subdf3$param == 'e.ssrd']), mean(subdf3$sd[subdf3$param == 'e.t2m']), mean(subdf3$sd[subdf3$param == 'e.tp']), mean(subdf3$sd[subdf3$param == 'e.utci']))

o.epi.mu3 <- c(mean(subdf3$mu[subdf3$param == 'o.absh']), mean(subdf3$mu[subdf3$param == 'o.ssrd']), mean(subdf3$mu[subdf3$param == 'o.t2m']), mean(subdf3$mu[subdf3$param == 'o.tp']), mean(subdf3$mu[subdf3$param == 'o.utci']))
o.epi.sd3 <- c(mean(subdf3$sd[subdf3$param == 'o.absh']), mean(subdf3$sd[subdf3$param == 'o.ssrd']), mean(subdf3$sd[subdf3$param == 'o.t2m']), mean(subdf3$sd[subdf3$param == 'o.tp']), mean(subdf3$sd[subdf3$param == 'o.utci']))


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
                        format.regtbl(e.ols.mu[4], e.ols.se[4]),
                        format.regtbl(e.ols.mu[5], e.ols.se[5]), rep(NA, 10)),
                  EPI.raw=c(format.regtbl(e.epi.mu2[1], e.epi.sd2[1]),
                            format.regtbl(e.epi.mu2[2], e.epi.sd2[2]),
                            format.regtbl(e.epi.mu2[3], e.epi.sd2[3]),
                            format.regtbl(e.epi.mu2[4], e.epi.sd2[4]),
                            format.regtbl(e.epi.mu2[5], e.epi.sd2[5]),
                            format.regtbl(o.epi.mu2[1], o.epi.sd2[1]),
                            format.regtbl(o.epi.mu2[2], o.epi.sd2[2]),
                            format.regtbl(o.epi.mu2[3], o.epi.sd2[3]),
                            format.regtbl(o.epi.mu2[4], o.epi.sd2[4]),
                            format.regtbl(o.epi.mu2[5], o.epi.sd2[5])),
                  EPI.country=c(format.regtbl(e.epi.mu3[1], e.epi.sd3[1]),
                                format.regtbl(e.epi.mu3[2], e.epi.sd3[2]),
                                format.regtbl(e.epi.mu3[3], e.epi.sd3[3]),
                                format.regtbl(e.epi.mu3[4], e.epi.sd3[4]),
                                format.regtbl(e.epi.mu3[5], e.epi.sd3[5]),
                                format.regtbl(o.epi.mu3[1], o.epi.sd3[1]),
                                format.regtbl(o.epi.mu3[2], o.epi.sd3[2]),
                                format.regtbl(o.epi.mu3[3], o.epi.sd3[3]),
                                format.regtbl(o.epi.mu3[4], o.epi.sd3[4]),
                                format.regtbl(o.epi.mu3[5], o.epi.sd3[5])),
                  EPI.hyper=c(format.regtbl(e.epi.mu[1], e.epi.sd[1]),
                              format.regtbl(e.epi.mu[2], e.epi.sd[2]),
                              format.regtbl(e.epi.mu[3], e.epi.sd[3]),
                              format.regtbl(e.epi.mu[4], e.epi.sd[4]),
                              format.regtbl(e.epi.mu[5], e.epi.sd[5]),
                              format.regtbl(o.epi.mu[1], o.epi.sd[1]),
                              format.regtbl(o.epi.mu[2], o.epi.sd[2]),
                              format.regtbl(o.epi.mu[3], o.epi.sd[3]),
                              format.regtbl(o.epi.mu[4], o.epi.sd[4]),
                              format.regtbl(o.epi.mu[5], o.epi.sd[5])))

names(tbl) <- c("OLS", "Epi (Raw)", "Epi (Country)", "Epi (Hyper)")
row.names(tbl) <- paste0(rep(c('e.', 'o.'), each=10), c(weather[1], " ", weather[2], "  ", weather[3], "   ", weather[4], "    ", weather[5], "     "))
xtable(tbl)

tbl1 <- tbl

tbl <- cbind(tbl1[, -2], tbl2[, -1:-2])
xtable(tbl)

