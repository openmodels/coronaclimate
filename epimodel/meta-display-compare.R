setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

library(dplyr)
library(ggplot2)
library(scales)
library(PBSmapping)

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

## Make raw histograms with mobile labelling

outfile <- "../../results/epimodel-meta-1217-all-nobs.csv"
mobilefile <- "../../results/epimodel-meta-1217-mobile-nobs.csv"
suffix <- "-1217-compare"

allrecorded <- read.csv(outfile)
mobrecorded <- read.csv(mobilefile)

results <- subset(allrecorded, group == "Raw")

results$mobile <- F
results$mobile[results$regid %in% unique(mobrecorded$regid)] <- T

results$param <- factor(results$param, paramorder)

## For plot, drop beyond the 99th
results$showit <- T
for (param in unique(results$param)) {
    limits <- quantile(results$mu[results$param == param], c(.025, .975), na.rm=T)
    results$showit[results$param == param] <- !is.na(results$mu[results$param == param]) & results$mu[results$param == param] > limits[1] & results$mu[results$param == param] < limits[2]
}
results$showit[results$param == 'mobility_slope' & !results$mobile] <- F

results$paramlabel <- as.character(results$param)
for (param in names(labelmap))
    results$paramlabel[results$param == param] <- labelmap[[param]]

results$paramlabel <- factor(results$paramlabel, levels=sapply(paramorder, function(param) labelmap[[param]]))

gp <- ggplot(subset(results, showit), aes(mu, fill=mobile)) +
    facet_wrap(~ paramlabel, scales='free') +
    geom_histogram() + xlab(NULL) + ylab(NULL) + theme_bw() +
    scale_fill_discrete(name="Mobility data") +
    theme(legend.justification=c(1,0), legend.position=c(.9,0))
ggsave(paste0("../../figures/raw-results", suffix, ".pdf"), width=12, height=7)


## Make violin plots

outfiles <- c("../../results/epimodel-meta-1217-all-nobs.csv", "../../results/epimodel-meta-1217-mobile-nobs.csv")
outlabels <- c("All", "Mobile only")
allrecorded <- data.frame()

for (cc in 1:length(outfiles)) {
    allrecsub <- read.csv(outfiles[cc])
    allrecsub$label <- outlabels[cc]
    allrecorded <- rbind(allrecorded, allrecsub)
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

ggplot(subset(allrecorded, Country == "" & paramgroup != "Drop"), aes(paste(paramlabel, label), mu)) +
    facet_wrap(~ paramgroup, ncol=1, scales="free") +
    coord_flip() +
    geom_violin(data=subset(allrecorded, Country != "" & group == "Raw" & Region == "" & paramgroup != "Drop"), colour="#83242480", scale="width") + # muted('red'), alpha=.5
    geom_violin(data=subset(allrecorded, Country != "" & group == "Combined" & Region == "" & paramgroup != "Drop"), fill=muted('blue'), alpha=.5, linetype='blank', scale="width") +
    geom_point() + geom_errorbar(aes(ymin=ci2.5, ymax=ci97.5)) +
    theme_bw() + ylab("Hyper-paramater value and 95% CI") + xlab(NULL)
