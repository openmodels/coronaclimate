setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

library(dplyr)
library(ggplot2)
library(scales)
library(PBSmapping)

outfile <- "../../results/epimodel-meta-0314-full-all-nobs-nodel.csv"
outfile.mob <- "../../results/epimodel-meta-0314-full-mobile-nobs-nodel.csv"
suffix <- "-0314-all"

do.old.figures <- F

allrecorded <- read.csv(outfile)

paramorder <- c('alpha', 'invgamma', 'invsigma', 'invkappa', 'invtheta',
                'mobility_slope', 'omega', 'portion_early',
                'deathrate', 'deathomegaplus', 'deathlearning',
                'logbeta', 'logomega', 'eein',
                'e.absh', 'e.r', 'e.t2m', 'e.tp', 'e.ssrd', 'e.utci',
                'o.absh', 'o.r', 'o.t2m', 'o.tp', 'o.ssrd', 'o.utci', 'error')
labelmap <- list('mobility_slope'="Mobility Adjustment",
                 'alpha'="Gradual Adjustment Rate",
                 'invsigma'="Incubation Period (days)",
                 'invgamma'="Infectious Period (days)",
                 'invkappa'="Pre-Testing Period (days)",
                 'invtheta'="Pre-Reporting Period (days)",
                 'portion_early'="Portion Detected Early",
                 'logbeta'="Log Transmission Rate",
                 'omega'="Recording Rate",
                 'deathrate'="Death Rate",
                 'deathomegaplus'="Extra Record of Deaths",
                 'deathlearning'="Death Learning Rate",
                 'portion_early'="Portion Reported Early",
                 'e.t2m'="Air Temperature Trans.",
                 'e.tp'="Total Precipitation Trans.",
                 'e.r'="Relative Humidity Trans.",
                 'e.absh'="Absolute Humidity Trans.",
                 'e.ssrd'="Solar Radiation Trans.",
                 'e.utci'="Thermal Discomfort Trans.",
                 'o.t2m'="Air Temperature Detect",
                 'o.tp'="Total Precipitation Detect",
                 'o.r'="Relative Humidity Detect",
                 'o.absh'="Absolute Humidity Detect",
                 'o.ssrd'="Solar Radiation Detect",
                 'o.utci'="Thermal Discomfort Detect",
                 'error'="Model Error", 'logomega'="Log Reporting Rate",
                 'eein'="Exposed Imports")

## Display raw results
results <- subset(allrecorded, group == "Raw")

results$param <- factor(results$param, paramorder)

## For plot, drop beyond the 99th
results$showit <- T
for (param in unique(results$param)) {
    limits <- quantile(results$mu[results$param == param], c(.025, .975), na.rm=T)
    results$showit[results$param == param] <- !is.na(results$mu[results$param == param]) & results$mu[results$param == param] > limits[1] & results$mu[results$param == param] < limits[2]
}

results$paramlabel <- as.character(results$param)
for (param in names(labelmap))
    results$paramlabel[results$param == param] <- labelmap[[param]]

results$paramlabel <- factor(results$paramlabel, levels=sapply(paramorder, function(param) labelmap[[param]]))

if (do.old.figures) {
    ggplot(subset(results, showit), aes(mu)) +
        facet_wrap(~ paramlabel, scales='free') +
        geom_histogram() + xlab(NULL) + ylab(NULL) + theme_bw()
    ggsave(paste0("../../figures/raw-results", suffix, ".pdf"), width=12, height=7)
}

if (do.old.figures) {
## Scatter plot
library(reshape2)
results2 <- dcast(subset(results, param != 'error' & showit), regid ~ param, value.var='mu')

# Correlation panel
panel.cor <- function(x, y){
    usr <- par("usr"); on.exit(par(usr))
    par(usr=c(0, 1, 0, 1))
    r <- round(cor(x, y, use='complete', method='spearman'), digits=2)
    txt <- as.character(r)
    text(0.5, 0.5, txt)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x, y, pch=19, cex=.01)
}

png("../../figures/scattermatrix.png", width=1000, height=1000)
pairs(results2[,-1], lower.panel=panel.cor, upper.panel=upper.panel)
dev.off()
}

## Show all results
for (pp in unique(allrecorded$param)) {
    ggplot(subset(allrecorded, param == pp & Region == "" & Locality == ""), aes(regid, mu, colour=group)) +
        coord_flip() +
        geom_point() + geom_errorbar(aes(ymin=ci25, ymax=ci75)) +
        theme_bw() + theme(axis.text.y=element_text(size=4)) + ylab(NULL) + xlab(NULL) +
        scale_colour_discrete(name=paste0(pp, "\nEstimates"))
    ggsave(paste0("../../figures/epimodel-param-", pp, suffix, ".png"), width=5, height=7)
}

allrecorded$paramlabel <- as.character(allrecorded$param)
for (param in names(labelmap))
    allrecorded$paramlabel[allrecorded$param == param] <- labelmap[[param]]

allrecorded$paramlabel <- factor(allrecorded$paramlabel, levels=rev(sapply(paramorder, function(param) labelmap[[param]])))

allrecorded$paramgroup <- "Drop"
allrecorded$paramgroup[allrecorded$param %in% c('invsigma', 'invkappa', 'invgamma', 'invtheta')] <- "Period Lengths"
allrecorded$paramgroup[allrecorded$param %in% c('e.t2m', 'e.tp', 'e.r', 'e.absh', 'e.ssrd', 'e.utci')] <- "Weather on Log Transmission"
allrecorded$paramgroup[allrecorded$param %in% c('o.t2m', 'o.tp', 'o.r', 'o.absh', 'o.ssrd', 'o.utci')] <- "Weather on Log Detection"
allrecorded$paramgroup[allrecorded$param %in% c('portion_early', 'omega', 'deathrate', 'deathomegaplus')] <- "Proportional Response"
allrecorded$paramgroup[allrecorded$param %in% c('mobility_slope', 'alpha')] <- "Behavioural Response"
allrecorded$paramgroup[allrecorded$param %in% c('logbeta', 'logomega', 'deathlearning')] <- "Baseline Log Rates"

for (paramgroup in unique(allrecorded$paramgroup)) {
    botlev <- quantile(allrecorded$mu[allrecorded$paramgroup == paramgroup], .01, na.rm=T)
    toplev <- quantile(allrecorded$mu[allrecorded$paramgroup == paramgroup], .99, na.rm=T)
    allrecorded$paramgroup[allrecorded$paramgroup == paramgroup & allrecorded$mu < botlev] <- "Drop"
    allrecorded$paramgroup[allrecorded$paramgroup == paramgroup & allrecorded$mu > toplev] <- "Drop"
}

allrecorded$paramgroup[grep("-\\d", allrecorded$Country)] <- "Drop"

## Also overlay mobility values
mobrecorded <- read.csv(outfile.mob)

mobrecorded$paramlabel <- as.character(mobrecorded$param)
for (param in names(labelmap))
    mobrecorded$paramlabel[mobrecorded$param == param] <- labelmap[[param]]

mobrecorded$paramlabel <- factor(mobrecorded$paramlabel, levels=rev(sapply(paramorder, function(param) labelmap[[param]])))

mobrecorded$paramgroup <- "Drop"
mobrecorded$paramgroup[mobrecorded$param %in% c('invsigma', 'invkappa', 'invgamma', 'invtheta')] <- "Period Lengths"
mobrecorded$paramgroup[mobrecorded$param %in% c('e.t2m', 'e.tp', 'e.r', 'e.absh', 'e.ssrd', 'e.utci')] <- "Weather on Log Transmission"
mobrecorded$paramgroup[mobrecorded$param %in% c('o.t2m', 'o.tp', 'o.r', 'o.absh', 'o.ssrd', 'o.utci')] <- "Weather on Log Detection"
mobrecorded$paramgroup[mobrecorded$param %in% c('portion_early', 'omega', 'deathrate', 'deathomegaplus')] <- "Proportional Response"
mobrecorded$paramgroup[mobrecorded$param %in% c('mobility_slope', 'alpha')] <- "Behavioural Response"
mobrecorded$paramgroup[mobrecorded$param %in% c('logbeta', 'logomega', 'deathlearning')] <- "Baseline Log Rates"

for (paramgroup in unique(mobrecorded$paramgroup)) {
    botlev <- quantile(mobrecorded$mu[mobrecorded$paramgroup == paramgroup], .01, na.rm=T)
    toplev <- quantile(mobrecorded$mu[mobrecorded$paramgroup == paramgroup], .99, na.rm=T)
    mobrecorded$paramgroup[mobrecorded$paramgroup == paramgroup & mobrecorded$mu < botlev] <- "Drop"
    mobrecorded$paramgroup[mobrecorded$paramgroup == paramgroup & mobrecorded$mu > toplev] <- "Drop"
}

mobrecorded$paramgroup[grep("-\\d", mobrecorded$Country)] <- "Drop"

allrecorded$paramgroup[allrecorded$param == "mobility_slope" & !(allrecorded$regid %in% unique(mobrecorded$regid))] <- "Drop"

global.mobrecorded <- subset(mobrecorded, Country == "" & paramgroup != "Drop")

allrecorded$paramgroup <- factor(allrecorded$paramgroup, c("Period Lengths", "Proportional Response", "Behavioural Response", "Baseline Log Rates", "Weather on Log Transmission", "Weather on Log Detection", "Drop"))
global.mobrecorded$paramgroup <- factor(global.mobrecorded$paramgroup, c("Period Lengths", "Proportional Response", "Behavioural Response", "Baseline Log Rates", "Weather on Log Transmission", "Weather on Log Detection", "Drop"))


ggplot(subset(allrecorded, Country == "" & paramgroup != "Drop"), aes(paramlabel, mu)) +
    facet_wrap(~ paramgroup, ncol=1, scales="free") +
    coord_flip() +
    geom_violin(data=subset(allrecorded, Country != "" & group == "Raw" & Region == "" & paramgroup != "Drop"), colour="#83242480", scale="width") + # muted('red'), alpha=.5
    geom_violin(data=subset(allrecorded, Country != "" & group == "Combined" & Region == "" & paramgroup != "Drop"), fill=muted('blue'), alpha=.5, linetype='blank', scale="width") +
        geom_point(data=global.mobrecorded, aes(colour='Mobility only')) + geom_errorbar(data=global.mobrecorded, aes(ymin=ci2.5, ymax=ci97.5, colour='Mobility only')) +
        geom_point(aes(colour='All observations')) + geom_errorbar(aes(ymin=ci2.5, ymax=ci97.5, colour='All observations')) +
        scale_colour_discrete(name=NULL) +
        theme_bw() + ylab("Hyper-paramater value and 95% CI") + xlab(NULL)
ggsave(paste0("../../figures/epimodel-violins", suffix, ".pdf"), width=8, height=7.5)

validgroups <- levels(allrecorded$paramgroup)[!(levels(allrecorded$paramgroup) %in% c("Drop", "Behavioural Response", "Baseline Log Rates"))]

ggplot(subset(allrecorded, Country == "" & paramgroup %in% validgroups), aes(paramlabel, mu)) +
    facet_wrap(~ paramgroup, ncol=1, scales="free") +
    coord_flip() +
    geom_violin(data=subset(allrecorded, Country != "" & group == "Raw" & Region == "" & paramgroup %in% validgroups), colour="#83242480", scale="width") + # muted('red'), alpha=.5
    geom_violin(data=subset(allrecorded, Country != "" & group == "Combined" & Region == "" & paramgroup %in% validgroups), fill=muted('blue'), alpha=.5, linetype='blank', scale="width") +
        geom_point(data=subset(global.mobrecorded, paramgroup %in% validgroups), aes(colour='Mobility only')) +
        geom_errorbar(data=subset(global.mobrecorded, paramgroup %in% validgroups), aes(ymin=ci2.5, ymax=ci97.5, colour='Mobility only')) +
        geom_point(aes(colour='All observations')) + geom_errorbar(aes(ymin=ci2.5, ymax=ci97.5, colour='All observations')) +
        scale_colour_discrete(name=NULL) +
        theme_bw() + ylab("Hyper-paramater value and 95% CI") + xlab(NULL)
ggsave(paste0("../../figures/epimodel-violins", suffix, "-present.pdf"), width=8, height=5.5)

## Look at variance across regions
library(xtable)

sumtbl <- subset(allrecorded, Country != "" & group == "Combined" & Region == "" & paramgroup != "Drop") %>% group_by(paramlabel) %>% summarize(mean.mu=mean(mu), sd.mu=sd(mu), mu.sd=mean(sd), ci.25=quantile(mu, .25), ci.75=quantile(mu, .75)) %>% left_join(subset(allrecorded, Country == "" & paramgroup != "Drop")[, c('paramlabel', 'mu')]) %>% left_join(subset(mobrecorded, Country == "" & paramgroup != "Drop")[, c('paramlabel', 'mu')], by='paramlabel')
sumtbl$sd.ratio <- sumtbl$sd.mu / sumtbl$mu.sd
print(xtable(sumtbl[nrow(sumtbl):1, c(1:2, 5:6, 9, 8, 7)], digits=3), include.rownames=F)


## Prepare to map

## Grab ALPHA.3 from df
df <- read.csv("../../cases/panel_all.csv")
df2 <- df[, c('Country', 'ALPHA.3')] %>% group_by(Country) %>% summarize(ALPHA.3=ALPHA.3[1])
allrecorded2 <- subset(allrecorded, Region == "" & Locality == "" & group == "Combined") %>% left_join(df2)

## Map to PIDs
shp <- importShapefile("../../shapefiles/gadm36_levels_simple/adm0.shp")
polydata <- attr(shp, 'PolyData')

allrecorded3 <- allrecorded2 %>% left_join(polydata, by=c('ALPHA.3'='GID_0'))

for (pp in unique(allrecorded$param)) {
    shp2 <- shp %>% left_join(allrecorded3[allrecorded3$param == pp, c('PID', 'mu')])

    ggplot(shp2, aes(X, Y, fill=mu, group=paste(PID, SID))) +
        geom_polygon() + scale_y_continuous(name=NULL, limits=c(-60, 85), expand=c(0, 0)) +
        scale_x_continuous(name=NULL, expand=c(0, 0)) + theme_bw() + scale_fill_continuous(name=pp)

    ggsave(paste0("../../figures/epimodel-param-map-", pp, ".png"), width=8, height=3)
}

if (do.old.figures) {
## Plot "total weather effect"

allsqssq <- subset(allrecorded3, param %in% c('e.absh', 'e.r', 'e.t2m', 'e.tp', 'e.ssrd', 'e.utci',
                                              'o.absh', 'o.r', 'o.t2m', 'o.tp', 'o.ssrd', 'o.utci')) %>% group_by(regid) %>% summarize(Country=Country[1], sssmu=sqrt(sum(mu^2)), e.absh=mu[param == 'e.absh'], e.t2m=mu[param == 'e.t2m'], e.tp=mu[param == 'e.tp'], e.ssrd=mu[param == 'e.ssrd'], e.utci=mu[param == 'e.utci'], o.absh=mu[param == 'o.absh'], o.t2m=mu[param == 'o.t2m'], o.tp=mu[param == 'o.tp'], o.ssrd=mu[param == 'o.ssrd'], o.utci=mu[param == 'o.utci'], PID=PID[1])
shp2 <- shp %>% left_join(allsqssq[, c('PID', 'sssmu')])

gp <- ggplot(shp2, aes(X, Y, fill=sssmu, group=paste(PID, SID))) +
    geom_polygon() + scale_y_continuous(name=NULL, limits=c(-60, 85), expand=c(0, 0)) +
    scale_x_continuous(name=NULL, expand=c(0, 0)) + theme_bw() + scale_fill_continuous(name="Weather effect", limits=c(0, max(shp2$sssmu)))
ggsave(paste0("../../figures/epimodel-weather-1217-all-sss.png"), width=8, height=3)

weatherscales <- apply(df[, weather], 2, sd)

df.var <- subset(df, Region == "" & Locality == "") %>% group_by(Country) %>% summarize(var.absh=var(absh / weatherscales[1]), var.ssrd=var(ssrd / weatherscales[2]), var.t2m=var(t2m / weatherscales[3]), var.tp=var(tp / weatherscales[3]), var.utci=var(utci / weatherscales[1]))

allsss2 <- allsqssq %>% left_join(df.var)
allsss2$totsss <- sqrt((allsss2$e.absh^2 + allsss2$o.absh^2) * allsss2$var.absh + (allsss2$e.t2m^2 + allsss2$o.t2m^2) * allsss2$var.t2m + (allsss2$e.tp^2 + allsss2$o.tp^2) * allsss2$var.tp + (allsss2$e.ssrd^2 + allsss2$o.ssrd^2) * allsss2$var.ssrd + (allsss2$e.utci^2 + allsss2$o.utci^2) * allsss2$var.utci)
shp2 <- shp %>% left_join(allsss2[, c('PID', 'totsss')])

gp <- ggplot(shp2, aes(X, Y, fill=totsss, group=paste(PID, SID))) +
    geom_polygon() + scale_y_continuous(name=NULL, limits=c(-60, 85), expand=c(0, 0)) +
    scale_x_continuous(name=NULL, expand=c(0, 0)) + theme_bw() + scale_fill_continuous(name="Weather variance", limits=c(0, max(shp2$totsss)))
ggsave(paste0("../../figures/epimodel-weather-1217-all-totsss.png"), width=8, height=3)
}

## Version 2: Estiamte the variation in logbeta and logomega from weather

library(lfe)

allparams <- subset(allrecorded3, param %in% c('logbeta', 'e.absh', 'e.r', 'e.t2m', 'e.tp', 'e.ssrd', 'e.utci',
                                               'logomega', 'o.absh', 'o.r', 'o.t2m', 'o.tp', 'o.ssrd', 'o.utci')) %>% group_by(regid) %>% summarize(Country=Country[1], logbeta=mean(mu[param == 'logbeta']), e.absh=mean(mu[param == 'e.absh']), e.t2m=mean(mu[param == 'e.t2m']), e.tp=mean(mu[param == 'e.tp']), e.ssrd=mean(mu[param == 'e.ssrd']), e.utci=mean(mu[param == 'e.utci']), logomega=mean(mu[param == 'logomega']), o.absh=mean(mu[param == 'o.absh']), o.t2m=mean(mu[param == 'o.t2m']), o.tp=mean(mu[param == 'o.tp']), o.ssrd=mean(mu[param == 'o.ssrd']), o.utci=mean(mu[param == 'o.utci']), PID=PID[1])

df$regid <- paste(df$Country, df$Region, df$Locality)

weather <- c('absh', 't2m', 'tp', 'ssrd', 'utci')
weatherscales <- apply(df[, weather], 2, sd)

plotdf <- data.frame()
for (regid in unique(allparams$regid)) {
    subdf <- df[df$regid == regid,]
    values <- demeanlist(subdf[, weather], list(factor(rep('all', nrow(subdf))))) / t(matrix(weatherscales, ncol=nrow(subdf), nrow=length(weather)))

    daily <- allparams$logbeta[allparams$regid == regid] + as.matrix(values) %*% t(allparams[allparams$regid == regid, paste0('e.', weather)])
    e.clfrac <- quantile(exp(daily), .025) / mean(exp(daily))
    e.chfrac <- quantile(exp(daily), .975) / mean(exp(daily))
    e.sdfrac <- sd(exp(daily)) / mean(exp(daily))

    omega <- exp(allparams$logomega[allparams$regid == regid])
    daily <- (omega / (1 + omega)) * exp(as.matrix(values) %*% t(allparams[allparams$regid == regid, paste0('o.', weather)]))
    o.clfrac <- quantile(daily, .025) / mean(daily)
    o.chfrac <- quantile(daily, .975) / mean(daily)
    o.sdfrac <- sd(daily) / mean(daily)

    plotdf <- rbind(plotdf, data.frame(PID=allparams$PID[allparams$regid == regid],
                                       e.sdfrac, o.sdfrac))
}

shp2 <- shp %>% left_join(plotdf)

gp <- ggplot(shp2, aes(X, Y, fill=e.sdfrac, group=paste(PID, SID))) +
    geom_polygon() + scale_y_continuous(name=NULL, limits=c(-60, 85), expand=c(0, 0)) +
    scale_x_continuous(name=NULL, expand=c(0, 0)) + theme_bw() + scale_fill_continuous(name="Beta SD\nChange", trans='log10', labels=scales::percent) +#limits=c(0, max(shp2$e.sdfrac))) +
    theme(legend.justification=c(0,0), legend.position=c(0.01,0.01))
ggsave(paste0("../../figures/epimodel-weather-1217-all-esd.png"), width=8, height=3)

gp <- ggplot(shp2, aes(X, Y, fill=o.sdfrac, group=paste(PID, SID))) +
    geom_polygon() + scale_y_continuous(name=NULL, limits=c(-60, 85), expand=c(0, 0)) +
    scale_x_continuous(name=NULL, expand=c(0, 0)) + theme_bw() + scale_fill_continuous(name="Omega SD\nChange", trans='log10', labels=scales::percent) +
    theme(legend.justification=c(0,0), legend.position=c(0.01,0.01))
ggsave(paste0("../../figures/epimodel-weather-1217-all-osd.png"), width=8, height=3)

mean(plotdf$e.sdfrac[!is.na(plotdf$PID)]) * 100
mean(plotdf$o.sdfrac[!is.na(plotdf$PID)]) * 100

quantile(plotdf$e.sdfrac[!is.na(plotdf$PID)], c(.025, .975)) * 100
quantile(plotdf$o.sdfrac[!is.na(plotdf$PID)], c(.025, .975)) * 100

## Plot the US

allrecorded2 <- subset(allrecorded, Country == "USA" & Locality == "" & group == "Combined")

shp <- importShapefile("../../shapefiles/gadm36_USA_shp/gadm36_USA_1_simple.shp")
polydata <- attr(shp, 'PolyData')

polydata2 <- polydata %>% left_join(data.frame(NAME_1=state.name, Region=state.abb))
allrecorded3 <- allrecorded2 %>% left_join(polydata2[, c('Region', 'PID')])

for (pp in unique(allrecorded$param)) {
    shp2 <- shp %>% left_join(allrecorded3[allrecorded3$param == pp, c('PID', 'mu')])

    ggplot(shp2, aes(X, Y, fill=mu, group=paste(PID, SID))) +
        geom_polygon() + scale_y_continuous(name=NULL, limits=c(25, 50)) +
        scale_x_continuous(name=NULL, limits=c(-125, -67)) + theme_bw() + scale_fill_continuous(name=pp)

    ggsave(paste0("../../figures/epimodel-param-usmap-", pp, ".png"), width=8, height=4)
}
