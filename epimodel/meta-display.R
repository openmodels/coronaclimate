setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

library(dplyr)
library(ggplot2)
library(scales)
library(PBSmapping)

outfile <- "../../results/epimodel-meta-0806-x5.csv"
allrecorded <- read.csv(outfile)

## Show all results
for (pp in unique(allrecorded$param)) {
    ggplot(subset(allrecorded, param == pp & Region == "" & Locality == ""), aes(regid, mu, colour=group)) +
        coord_flip() +
        geom_point() + geom_errorbar(aes(ymin=ci25, ymax=ci75)) +
        theme_bw() + theme(axis.text.y=element_text(size=4)) + ylab(NULL) + xlab(NULL) +
        scale_colour_discrete(name=paste0(pp, "\nEstimates"))
    ggsave(paste0("../../figures/epimodel-param-", pp, ".png"), width=5, height=7)
}

allrecorded$paramgroup <- "Drop"
allrecorded$paramgroup[allrecorded$param %in% c('invsigma', 'invkappa', 'invgamma')] <- "Period Lengths"
allrecorded$paramgroup[allrecorded$param %in% c('e.tp', 'e.r', 'e.absh')] <- "Weather Response"
allrecorded$paramgroup[allrecorded$param %in% c('omega', 'deathrate', 'deathomegaplus')] <- "Proportional Response"
allrecorded$paramgroup[allrecorded$param %in% c('mobility_slope', 'portion_early', 'alpha')] <- "Behavioural Response"

for (paramgroup in unique(allrecorded$paramgroup)) {
    botlev <- quantile(allrecorded$mu[allrecorded$paramgroup == paramgroup], .01)
    toplev <- quantile(allrecorded$mu[allrecorded$paramgroup == paramgroup], .99)
    allrecorded$paramgroup[allrecorded$paramgroup == paramgroup & allrecorded$mu < botlev] <- "Drop"
    allrecorded$paramgroup[allrecorded$paramgroup == paramgroup & allrecorded$mu > toplev] <- "Drop"
}

ggplot(subset(allrecorded, Country == "" & paramgroup != "Drop"), aes(param, mu)) +
    facet_wrap(~ paramgroup, ncol=1, scales="free") +
    coord_flip() +
    geom_violin(data=subset(allrecorded, Country != "" & Region == "" & paramgroup != "Drop"), fill=muted('blue'), alpha=.5, linetype='blank', scale="width") +
    geom_point() + geom_errorbar(aes(ymin=ci2.5, ymax=ci97.5)) +
    theme_bw() + ylab("Hyper-paramater value and 95% CI") + xlab(NULL)

## Prepare to map

## Grab ALPHA.3 from df
load("../../cases/panel-prepped_MLI.RData")
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


