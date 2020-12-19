setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

mainmodel <- '1201'

weather <- c('absh', 'ssrd', 't2m', 'tp', 'utci')
e.ols.mu <- c(0.03344, -0.00332, -0.04791, -0.00190, -0.00684)
e.ols.se <- c(0.01493, 0.00423, 0.01812, 0.00135, 0.00632)

alltbls <- list()
allgdfs <- list()

for (filepath in c('../../results/epimodel-meta-1030-all-pop.csv', '../../results/epimodel-meta-1030-mobile-pop.csv', paste0('../../results/epimodel-meta-', mainmodel, '-all-nobs.csv'), paste0('../../results/epimodel-meta-', mainmodel, '-mobile-nobs.csv'))) {

df <- read.csv(filepath)
subdf <- subset(df, Country == "" & Region == "")

e.epi.mu <- c(subdf$mu[subdf$param == 'e.absh'], subdf$mu[subdf$param == 'e.ssrd'], subdf$mu[subdf$param == 'e.t2m'], subdf$mu[subdf$param == 'e.tp'], subdf$mu[subdf$param == 'e.utci'])
e.epi.sd <- c(subdf$sd[subdf$param == 'e.absh'], subdf$sd[subdf$param == 'e.ssrd'], subdf$sd[subdf$param == 'e.t2m'], subdf$sd[subdf$param == 'e.tp'], subdf$sd[subdf$param == 'e.utci'])

o.epi.mu <- c(subdf$mu[subdf$param == 'o.absh'], subdf$mu[subdf$param == 'o.ssrd'], subdf$mu[subdf$param == 'o.t2m'], subdf$mu[subdf$param == 'o.tp'], subdf$mu[subdf$param == 'o.utci'])
o.epi.sd <- c(subdf$sd[subdf$param == 'o.absh'], subdf$sd[subdf$param == 'o.ssrd'], subdf$sd[subdf$param == 'o.t2m'], subdf$sd[subdf$param == 'o.tp'], subdf$sd[subdf$param == 'o.utci'])

allgdfs[[filepath]] <- data.frame(weather=rep(weather, 2), mu=c(e.epi.mu, o.epi.mu), sd=c(e.epi.sd, o.epi.sd), channel=rep(c('Transmission', 'Detection'), each=length(weather)))

subdf2 <- subset(df, Country != "" & Region == "" & group == "Raw")
e.epi.mu2 <- c(mean(subdf2$mu[subdf2$param == 'e.absh']), mean(subdf2$mu[subdf2$param == 'e.ssrd']), mean(subdf2$mu[subdf2$param == 'e.t2m']), mean(subdf2$mu[subdf2$param == 'e.tp']), mean(subdf2$mu[subdf2$param == 'e.utci']))
e.epi.sd2 <- c(mean(subdf2$sd[subdf2$param == 'e.absh']), mean(subdf2$sd[subdf2$param == 'e.ssrd']), mean(subdf2$sd[subdf2$param == 'e.t2m']), mean(subdf2$sd[subdf2$param == 'e.tp']), mean(subdf2$sd[subdf2$param == 'e.utci']))

o.epi.mu2 <- c(mean(subdf2$mu[subdf2$param == 'o.absh']), mean(subdf2$mu[subdf2$param == 'o.ssrd']), mean(subdf2$mu[subdf2$param == 'o.t2m'], na.rm=T), mean(subdf2$mu[subdf2$param == 'o.tp']), mean(subdf2$mu[subdf2$param == 'o.utci']))
o.epi.sd2 <- c(mean(subdf2$sd[subdf2$param == 'o.absh']), mean(subdf2$sd[subdf2$param == 'o.ssrd']), mean(subdf2$sd[subdf2$param == 'o.t2m'], na.rm=T), mean(subdf2$sd[subdf2$param == 'o.tp']), mean(subdf2$sd[subdf2$param == 'o.utci']))

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

alltbls[[filepath]] <- tbl
}

## Produce table
tbl.row.names <- sapply(0:19, function(ii) ifelse(ii < 10, ifelse(ii %% 2 == 0, paste("Transmission", weather[1 + ii / 2]), ""), ifelse(ii %% 2 == 0, paste("Detection", weather[1 + (ii - 10) / 2]), "")))
tbl <- cbind(tbl.row.names, alltbls[[paste0('../../results/epimodel-meta-', mainmodel, '-mobile-nobs.csv')]][, -2],
             alltbls[[paste0('../../results/epimodel-meta-', mainmodel, '-all-nobs.csv')]][, -1:-2])
names(tbl) <- c("", "OLS", "Bayes (Country)", "Bayes (Hyper)", "Bayes (Country)", "Bayes (Hyper)")

print(xtable(tbl), include.rownames=F)

## Produce bars

gdf <- rbind(data.frame(weather=rep(weather, 2), mu=rep(e.ols.mu, 2), sd=rep(e.ols.se, 2), channel="OLS", panel=rep(c("Mobility-Only", "All Observations"), each=5)),
             cbind(panel="Mobility-Only", allgdfs[[paste0('../../results/epimodel-meta-', mainmodel, '-mobile-nobs.csv')]]),
             cbind(panel="All Observations", allgdfs[[paste0('../../results/epimodel-meta-', mainmodel, '-all-nobs.csv')]]))

library(ggplot2)

gdf$panel <- factor(gdf$panel, levels=c('Mobility-Only', 'All Observations'))
gdf$channel <- factor(gdf$channel, levels=c('OLS', 'Transmission', 'Detection'))

ggplot(gdf, aes(weather, mu, fill=channel)) +
    facet_grid(. ~ panel) +
    geom_bar(stat='identity', position='dodge') + geom_errorbar(aes(ymin=mu - 1.96*sd, ymax=mu + 1.96*sd), position=position_dodge(width=.875), width=.5) + theme_bw() + xlab(NULL) + ylab(NULL) +
    scale_fill_discrete(name="Channel:")
