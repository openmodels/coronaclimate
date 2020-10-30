setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

library(dplyr)
library(ggplot2)
library(scales)
library(PBSmapping)

outfile <- "../../results/epimodel-meta-1018-pop.csv"
allrecorded <- read.csv(outfile)

paramorder <- c('alpha', 'invgamma', 'invsigma', 'invkappa', 'invtheta',
                'mobility_slope', 'omega', 'portion_early',
                'deathrate', 'deathomegaplus', 'error',
                'logbeta', 'logomega', 'eein',
                'e.absh', 'e.r', 'e.t2m', 'e.tp', 'e.ssrd', 'e.utci',
                'o.absh', 'o.r', 'o.t2m', 'o.tp', 'o.ssrd', 'o.utci')
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

ggplot(subset(results, param != 'error' & showit), aes(mu)) +
    facet_wrap(~ paramlabel, scales='free') +
    geom_histogram() + xlab(NULL) + ylab(NULL) + theme_bw()

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

pairs(results2[,-1], lower.panel=panel.cor, upper.panel=upper.panel)

## Show all results
for (pp in unique(allrecorded$param)) {
    ggplot(subset(allrecorded, param == pp & Region == "" & Locality == ""), aes(regid, mu, colour=group)) +
        coord_flip() +
        geom_point() + geom_errorbar(aes(ymin=ci25, ymax=ci75)) +
        theme_bw() + theme(axis.text.y=element_text(size=4)) + ylab(NULL) + xlab(NULL) +
        scale_colour_discrete(name=paste0(pp, "\nEstimates"))
    ggsave(paste0("../../figures/epimodel-param-", pp, ".png"), width=5, height=7)
}

allrecorded$paramlabel <- as.character(allrecorded$param)
for (param in names(labelmap))
    allrecorded$paramlabel[allrecorded$param == param] <- labelmap[[param]]

allrecorded$paramlabel <- factor(allrecorded$paramlabel, levels=rev(sapply(paramorder, function(param) labelmap[[param]])))

allrecorded$paramgroup <- "Drop"
allrecorded$paramgroup[allrecorded$param %in% c('invsigma', 'invkappa', 'invgamma', 'invtheta')] <- "Period Lengths"
allrecorded$paramgroup[allrecorded$param %in% c('e.t2m', 'e.tp', 'e.r', 'e.absh', 'e.ssrd', 'e.utci')] <- "Weather on Transmission"
allrecorded$paramgroup[allrecorded$param %in% c('o.t2m', 'o.tp', 'o.r', 'o.absh', 'o.ssrd', 'o.utci')] <- "Weather on Detection"
allrecorded$paramgroup[allrecorded$param %in% c('portion_early', 'omega', 'deathrate', 'deathomegaplus')] <- "Proportional Response"
allrecorded$paramgroup[allrecorded$param %in% c('mobility_slope', 'alpha')] <- "Behavioural Response"

for (paramgroup in unique(allrecorded$paramgroup)) {
    botlev <- quantile(allrecorded$mu[allrecorded$paramgroup == paramgroup], .01, na.rm=T)
    toplev <- quantile(allrecorded$mu[allrecorded$paramgroup == paramgroup], .99, na.rm=T)
    allrecorded$paramgroup[allrecorded$paramgroup == paramgroup & allrecorded$mu < botlev] <- "Drop"
    allrecorded$paramgroup[allrecorded$paramgroup == paramgroup & allrecorded$mu > toplev] <- "Drop"
}

allrecorded$paramgroup[grep("-\\d", allrecorded$Country)] <- "Drop"

ggplot(subset(allrecorded, Country == "" & paramgroup != "Drop"), aes(paramlabel, mu)) +
    facet_wrap(~ paramgroup, ncol=1, scales="free") +
    coord_flip() +
    geom_violin(data=subset(allrecorded, Country != "" & group == "Raw" & Region == "" & paramgroup != "Drop"), colour="#83242480", scale="width") + # muted('red'), alpha=.5
    geom_violin(data=subset(allrecorded, Country != "" & group == "Combined" & Region == "" & paramgroup != "Drop"), fill=muted('blue'), alpha=.5, linetype='blank', scale="width") +
    geom_point() + geom_errorbar(aes(ymin=ci2.5, ymax=ci97.5)) +
    theme_bw() + ylab("Hyper-paramater value and 95% CI") + xlab(NULL)

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
