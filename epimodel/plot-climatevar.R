setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

library(dplyr)
library(reshape2)

dfout <- read.csv("../../results/epimodel-0105noprior-nodel.csv")
dfout$rank <- NA
for (param in unique(dfout$param))
    dfout$rank[dfout$param == param] <- rank(dfout$mu[dfout$param == param])

dfin <- read.csv("../../cases/panel_all.csv")
dfin$regid <- paste(dfin$Country, dfin$Region, dfin$Locality)

dfin2 <- dfin %>% group_by(regid) %>% summarize(ALPHA.3=ALPHA.3[1], Country=Country[1], Region=Region[1], Locality=Locality[1], area=mean(population / population_density), population_density=mean(population_density), nobs=length(population), lowest_level=min(lowest_level), implausible=max(implausible), absh=mean(absh), absh2=mean(absh2), ssrd=mean(ssrd), ssrd2=mean(ssrd2), t2m=mean(t2m), t2m2=mean(t2m2), tp=mean(tp), tp2=mean(tp2), utci=mean(utci), utci2=mean(utci2))

df <- dfout %>% left_join(dfin2)
df$logarea <- log(df$area)
df$logpopden <- log(df$population_density)

dfgov <- subset(read.csv("../../policy/governance_data/imputed_2019_governance_indicators.csv"), type == "Estimate")
dfgov2 <- dfgov %>% left_join(read.csv("../../policy/governance_data/meta_countries.csv"), by=c('code'='id'))

load("../../governance/WVS_Cross-National_Wave_7_R_v1_5.rdata")

unique(df$Country[df$ALPHA.3 == ""])
df$ALPHA.3[df$ALPHA.3 == "" & df$Country == "Australia"] <- "AUS"
df$ALPHA.3[df$ALPHA.3 == "" & df$Country == "Canada"] <- "CAN"
df2 <- df %>% left_join(dfgov2, by=c('ALPHA.3'='iso3'))

get.pattern <- function(param, subdf) {
    ## In levels
    mod <- lm(mu ~ logarea + logpopden + lowest_level + implausible + absh + absh2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 + utci + utci2 + incomeLevel + geo.lat + landlocked, data=subdf)
    rsqr <- summary(mod)$r.squared
    areave <- sum(anova(mod)[1, 2]) / sum(anova(mod)[, 2])
    popve <- sum(anova(mod)[2, 2]) / sum(anova(mod)[, 2])
    adminve1 <- sum(anova(mod)[3, 2]) / sum(anova(mod)[, 2])
    adminve2 <- sum(anova(mod)[4, 2]) / sum(anova(mod)[, 2])
    climve <- sum(anova(mod)[5:14, 2]) / sum(anova(mod)[, 2])
    incve <- sum(anova(mod)[15, 2]) / sum(anova(mod)[, 2])
    geove <- sum(anova(mod)[16:17, 2]) / sum(anova(mod)[, 2])

    mod <- lm(mu ~ logarea + logpopden + incomeLevel + lowest_level + implausible + absh + absh2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 + utci + utci2 + geo.lat + landlocked + reform + local + power_govern + hope_young + coordination + confidence + respect_law, data=subdf)
    rsqr.gov <- summary(mod)$r.squared

    mod <- lm(mu ~ factor(Country) + logarea + logpopden + lowest_level + implausible + absh + absh2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 + utci + utci2, data=subdf)
    rsqr.fe <- summary(mod)$r.squared

    row1 <- data.frame(param, values='levels', rsqr, areave, popve, incve, adminve1, adminve2, climve, geove, govve=rsqr.gov - rsqr, countryve=rsqr.fe - rsqr.gov)

    ## In ranks
    mod <- lm(rank ~ logarea + logpopden + lowest_level + implausible + absh + absh2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 + utci + utci2 + incomeLevel + geo.lat + landlocked, data=subdf)
    rsqr <- summary(mod)$r.squared
    areave <- sum(anova(mod)[1, 2]) / sum(anova(mod)[, 2])
    popve <- sum(anova(mod)[2, 2]) / sum(anova(mod)[, 2])
    adminve1 <- sum(anova(mod)[3, 2]) / sum(anova(mod)[, 2])
    adminve2 <- sum(anova(mod)[4, 2]) / sum(anova(mod)[, 2])
    climve <- sum(anova(mod)[5:14, 2]) / sum(anova(mod)[, 2])
    incve <- sum(anova(mod)[15, 2]) / sum(anova(mod)[, 2])
    geove <- sum(anova(mod)[16:17, 2]) / sum(anova(mod)[, 2])

    mod <- lm(rank ~ logarea + logpopden + incomeLevel + lowest_level + implausible + absh + absh2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 + utci + utci2 + geo.lat + landlocked + reform + local + power_govern + hope_young + coordination + confidence + respect_law, data=subdf)
    rsqr.gov <- summary(mod)$r.squared

    mod <- lm(rank ~ factor(Country) + logarea + logpopden + lowest_level + implausible + absh + absh2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 + utci + utci2, data=subdf)
    rsqr.fe <- summary(mod)$r.squared

    row2 <- data.frame(param, values='ranks', rsqr, areave, popve, incve, adminve1, adminve2, climve, geove, govve=rsqr.gov - rsqr, countryve=rsqr.fe - rsqr.gov)

    rbind(row1, row2)
}

patterns <- data.frame()
for (param in unique(df2$param)) {
    subdf <- df2[df2$param == param,]

    patterns <- rbind(patterns, get.pattern(param, subdf))
}

subdf <- subset(df2, param == 'alpha') %>% left_join(df2 %>% group_by(regid) %>% summarize(murhat=mean(rhat, na.rm=T)))
subdf$rank <- subdf$murhat
subdf$rank <- rank(subdf$murhat)
patterns <- rbind(patterns, get.pattern('rhat', subdf))

library(ggplot2)

patterns2 <- subset(melt(patterns, id.vars=c('param', 'values')), variable != 'rsqr')
patterns2$label <- NA
patterns2$label[patterns2$variable == 'areave'] <- "Area"
patterns2$label[patterns2$variable == 'popve'] <- "Pop. Density"
patterns2$label[patterns2$variable == 'adminve1'] <- "Lowest Level Flag"
patterns2$label[patterns2$variable == 'adminve2'] <- "Implausible Flag"
patterns2$label[patterns2$variable == 'climve'] <- "Climate (quadratics)"
patterns2$label[patterns2$variable == 'incve'] <- "Income level (Country)"
patterns2$label[patterns2$variable == 'geove'] <- "Geography (Country)"
patterns2$label[patterns2$variable == 'govve'] <- "Governance (Country)"
patterns2$label[patterns2$variable == 'countryve'] <- "Other country-specific"
patterns2$label <- factor(patterns2$label, levels=rev(c("Area", "Pop. Density", "Lowest Level Flag", "Implausible Flag", "Climate (quadratics)", "Income level (Country)", "Geography (Country)", "Governance (Country)", "Other country-specific")))

source("plotlib.R")

patterns2$param.label <- NA
for (ii in 1:nrow(patterns2))
    patterns2$param.label[ii] <- labelmap[[patterns2$param[ii]]]
patterns2$param.label <- factor(patterns2$param.label, levels=rev(sapply(paramorder, function(param) labelmap[[param]])))

ggplot(patterns2, aes(param.label, value, fill=label)) +
    facet_wrap(~ values) +
    coord_flip() + xlab(NULL) + ylab("Explained variance") +
    geom_col() + scale_y_continuous(labels=scales::percent) +
    scale_fill_discrete(name="Source of\nvariance:", breaks=c("Area", "Pop. Density", "Lowest Level Flag", "Implausible Flag", "Climate (quadratics)", "Income level (Country)", "Geography (Country)", "Governance (Country)", "Other country-specific")) + theme_bw()

patterns %>% group_by(values) %>% summarize(rsqr=mean(rsqr),
                                            areave=mean(areave),
                                            popve=mean(popve),
                                            incve=mean(incve),
                                            adminve1=mean(adminve1),
                                            adminve2=mean(adminve2),
                                            climve=mean(climve),
                                            geove=mean(geove),
                                            govve=mean(govve),
                                            countryve=mean(countryve))

## tree <- rpart(mu ~ population + lowest_level + implausible + absh + absh2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 + utci + utci2, data=subset(df,


## library(rpart)
## library(rpart.plot)



