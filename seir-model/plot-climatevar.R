setwd("~/research/coronavirus/code/seir-model")

library(dplyr)
library(reshape2)

dfout <- read.csv("../../results/epimodel-0314-noprior-nodel.csv")
dfout$rank <- NA
for (param in unique(dfout$param))
    dfout$rank[dfout$param == param] <- rank(dfout$mu[dfout$param == param])

dfin <- read.csv("../../cases/panel_all.csv")
dfin$regid <- paste(dfin$Country, dfin$Region, dfin$Locality)

dfin2 <- dfin %>% group_by(regid) %>% summarize(ALPHA.3=ALPHA.3[1], Country=Country[1], Region=Region[1], Locality=Locality[1], area=mean(population / population_density), population_density=mean(population_density), nobs=length(population), lowest_level=min(lowest_level), implausible=max(implausible), absh=mean(absh), absh2=mean(absh2), ssrd=mean(ssrd), ssrd2=mean(ssrd2), t2m=mean(t2m), t2m2=mean(t2m2), tp=mean(tp), tp2=mean(tp2), utci=mean(utci), utci2=mean(utci2))

df <- dfout %>% left_join(dfin2)
df$logarea <- log(df$area)
df$logpopden <- log(df$population_density)

dfgov <- subset(read.csv("~/Dropbox/Coronavirus and Climate/policy/governance_data/imputed_2019_governance_indicators.csv"), type == "Estimate")
dfgov2 <- dfgov %>% left_join(read.csv("~/Dropbox/Coronavirus and Climate/policy/governance_data/meta_countries.csv"), by=c('code'='id'))

dfge <- read.csv("~/Dropbox/Coronavirus and Climate/governance/wgi_2019.csv")
dfge$ge <- as.numeric(dfge$ge)

df$ALPHA.3[df$ALPHA.3 == "" & df$Country == "Australia"] <- "AUS"
df$ALPHA.3[df$ALPHA.3 == "" & df$Country == "Canada"] <- "CAN"
unique(df$Country[df$ALPHA.3 == ""])

df2 <- df %>% filter(lowest_level == 1) %>% left_join(dfgov2, by=c('ALPHA.3'='iso3')) %>% left_join(dfge, by='ALPHA.3', suffix=c('', '.ge'))

get.pattern <- function(param, subdf) {
    ## In levels
    mod <- lm(mu ~ logarea + logpopden + incomeLevel, data=subdf)
    rsqr.precl <- summary(mod)$r.squared
    areave <- sum(anova(mod)[1, 2]) / sum(anova(mod)[, 2])
    popve <- sum(anova(mod)[2, 2]) / sum(anova(mod)[, 2])
    ##adminve <- sum(anova(mod)[3, 2]) / sum(anova(mod)[, 2])
    incve <- sum(anova(mod)[3, 2]) / sum(anova(mod)[, 2])
    ##geove <- sum(anova(mod)[14:15, 2]) / sum(anova(mod)[, 2])
    rss.precl <- sum(mod$residuals^2) # 6 coeff

    mod <- lm(mu ~ logarea + logpopden + absh + absh2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 + utci + utci2 + incomeLevel, data=subdf)
    rsqr <- summary(mod)$r.squared
    climve <- sum(anova(mod)[3:12, 2]) / sum(anova(mod)[, 2])
    rss <- sum(mod$residuals^2) # 16 coeff

    mod <- lm(mu ~ logarea + logpopden + incomeLevel + absh + absh2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 + utci + utci2 + ge, data=subdf)
    rsqr.ge <- summary(mod)$r.squared
    rss.ge <- sum(mod$residuals^2) # 17 coeff

    mod <- lm(mu ~ logarea + logpopden + incomeLevel + absh + absh2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 + utci + utci2 + ge + reform + local + power_govern + hope_young + coordination + confidence + respect_law, data=subdf)
    rsqr.gov <- summary(mod)$r.squared
    rss.gov <- sum(mod$residuals^2) # 24 coeff

    mod <- lm(mu ~ factor(Country) + logarea + logpopden + absh + absh2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 + utci + utci2, data=subdf)
    rsqr.fe <- summary(mod)$r.squared
    rss.fe <- sum(mod$residuals^2) # 167 coeff

    row1 <- data.frame(param, values='levels', rsqr, areave, popve, incve, climve, geve=rsqr.ge - rsqr, govve=rsqr.gov - rsqr.ge, countryve=rsqr.fe - rsqr.gov, rss.precl, rss, rss.ge, rss.gov, rss.fe)

    ## In ranks
    mod <- lm(rank ~ logarea + logpopden + incomeLevel, data=subdf)
    rsqr.precl <- summary(mod)$r.squared
    areave <- sum(anova(mod)[1, 2]) / sum(anova(mod)[, 2])
    popve <- sum(anova(mod)[2, 2]) / sum(anova(mod)[, 2])
    ##adminve <- sum(anova(mod)[3, 2]) / sum(anova(mod)[, 2])
    incve <- sum(anova(mod)[3, 2]) / sum(anova(mod)[, 2])
    ##geove <- sum(anova(mod)[14:15, 2]) / sum(anova(mod)[, 2])
    rss.precl <- sum(mod$residuals^2) # 6 coeff

    mod <- lm(rank ~ logarea + logpopden + absh + absh2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 + utci + utci2 + incomeLevel, data=subdf)
    rsqr <- summary(mod)$r.squared
    climve <- sum(anova(mod)[3:12, 2]) / sum(anova(mod)[, 2])
    rss <- sum(mod$residuals^2) # 16 coeff

    mod <- lm(rank ~ logarea + logpopden + incomeLevel + absh + absh2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 + utci + utci2 + ge, data=subdf)
    rsqr.ge <- summary(mod)$r.squared
    rss.ge <- sum(mod$residuals^2) # 17 coeff

    mod <- lm(rank ~ logarea + logpopden + incomeLevel + absh + absh2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 + utci + utci2 + ge + reform + local + power_govern + hope_young + coordination + confidence + respect_law, data=subdf)
    rsqr.gov <- summary(mod)$r.squared
    rss.gov <- sum(mod$residuals^2) # 24 coeff

    mod <- lm(rank ~ factor(Country) + logarea + logpopden + absh + absh2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 + utci + utci2, data=subdf)
    rsqr.fe <- summary(mod)$r.squared
    rss.fe <- sum(mod$residuals^2) # 167 coeff

    row2 <- data.frame(param, values='ranks', rsqr, areave, popve, incve, climve, geve=rsqr.ge - rsqr, govve=rsqr.gov - rsqr.ge, countryve=rsqr.fe - rsqr.gov, rss.precl, rss, rss.ge, rss.gov, rss.fe)

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

patterns2 <- subset(melt(patterns, id.vars=c('param', 'values')), !(variable %in% c('rsqr', 'rss.precl', 'rss', 'rss.ge', 'rss.gov', 'rss.fe')))
patterns2$label <- NA
patterns2$label[patterns2$variable == 'areave'] <- "Area"
patterns2$label[patterns2$variable == 'popve'] <- "Pop. Density"
patterns2$label[patterns2$variable == 'adminve'] <- "Implausible Flag"
patterns2$label[patterns2$variable == 'climve'] <- "Climate (quadratics)"
patterns2$label[patterns2$variable == 'incve'] <- "Income level (Country)"
patterns2$label[patterns2$variable == 'geove'] <- "Geography (Country)"
patterns2$label[patterns2$variable == 'geve'] <- "Governance (Country)" #"Effectiveness (Country)"
patterns2$label[patterns2$variable == 'govve'] <- "Governance (Country)" #"Other country-specific"
patterns2$label[patterns2$variable == 'countryve'] <- "Other country-specific"
patterns2$label <- factor(patterns2$label, levels=rev(c("Area", "Pop. Density", "Lowest Level Flag", "Implausible Flag", "Climate (quadratics)", "Income level (Country)", "Geography (Country)", "Effectiveness (Country)", "Governance (Country)", "Other country-specific")))

source("plotlib.R")

patterns2$param.label <- get.param.labels(patterns2$param)

ggplot(patterns2, aes(param.label, value, fill=label)) +
    facet_wrap(~ values) +
    coord_flip() + xlab(NULL) + ylab("Explained variance") +
    geom_col() + scale_y_continuous(labels=scales::percent) +
    scale_fill_discrete(name="Source of variance:", breaks=c("Area", "Pop. Density", "Lowest Level Flag", "Implausible Flag", "Climate (quadratics)", "Income level (Country)", "Geography (Country)", "Effectiveness (Country)", "Governance (Country)", "Other country-specific")) + theme_bw()
ggsave("climatevar-0314.pdf", width=10, height=5)

patterns %>% group_by(values) %>% summarize(rsqr=mean(rsqr),
                                            areave=mean(areave),
                                            popve=mean(popve),
                                            incve=mean(incve),
                                            ##adminve=mean(adminve),
                                            climve=mean(climve),
                                            ##geove=mean(geove),
                                            geve=mean(geve),
                                            govve=mean(govve),
                                            countryve=mean(countryve))

## Calculate f-tests

get.pattern.full <- function(param, subdf) {
    mod <- lm(mu ~ factor(Country) + logarea + logpopden, data=subdf)
    rss.fenocl <- sum(mod$residuals^2) # 157 coeff

    mod <- lm(mu ~ factor(Country) + logarea + logpopden + absh + absh2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 + utci + utci2, data=subdf)
    rss.fewicl <- sum(mod$residuals^2) # 167 coeff

    mod <- lm(mu ~ logarea + logpopden + incomeLevel + absh + absh2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 + utci + utci2, data=subdf)
    rss.nogov <- sum(mod$residuals^2) # 14 coeff

    mod <- lm(mu ~ logarea + logpopden + incomeLevel + absh + absh2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 + utci + utci2 + ge + reform + local + power_govern + hope_young + coordination + confidence + respect_law, data=subdf)
    rss.wigov <- sum(mod$residuals^2) # 22 coeff

    data.frame(param, rss.fenocl, rss.fewicl, rss.nogov, rss.wigov, rss.error=sum(subdf$sd^2))
}

patterns <- data.frame()
for (param in unique(df2$param)) {
    subdf <- df2[df2$param == param,]
    patterns <- rbind(patterns, get.pattern.full(param, subdf))
}

paramsets <- list('trans'=c('invsigma', 'invgamma', 'mobility_slope', 'logbeta'),
                  'detect'=c('invkappa', 'invtheta', 'logomega'),
                  'deaths'=c('deathrate', 'deathlearning', 'deathomegaplus'),
                  'trans.weather'=c("e.t2m", "e.tp", "e.ssrd", "e.utci"),
                  'detect.weather'=c("o.t2m", "o.tp", "o.ssrd", "o.utci"))

fstats <- data.frame()
for (group in names(paramsets)) {
    ## Channel: Climate
    fstat.clim <- mean(((patterns$rss.fenocl - patterns$rss.fewicl)[patterns$param %in% paramsets[[group]]] / 10) / ((patterns$rss.fewicl + patterns$rss.error)[patterns$param %in% paramsets[[group]]] / (sum(df2$param == 'alpha') - 167)))
    ## Channel: Governance
    fstat.gov <- mean(((patterns$rss.nogov - patterns$rss.wigov)[patterns$param %in% paramsets[[group]]] / 8) / ((patterns$rss.wigov + patterns$rss.error)[patterns$param %in% paramsets[[group]]] / (sum(df2$param == 'alpha') - 22)))
    fstats <- rbind(fstats, data.frame(group, fstat.clim, fstat.gov))
}

fstats$pvals.clim <- (1 - pf(fstats$fstat.clim, 10, sum(df2$param == 'alpha') - 167))
fstats$pvals.gov <- (1 - pf(fstats$fstat.gov, 8, sum(df2$param == 'alpha') - 22))


## tree <- rpart(mu ~ population + absh + absh2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 + utci + utci2, data=subset(df,


## library(rpart)
## library(rpart.plot)



