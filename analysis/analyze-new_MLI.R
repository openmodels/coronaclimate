library(dplyr)
library(lfe)
library(ggplot2)
library(stargazer)
setwd("~/Dropbox/Coronavirus and Climate")

SWITCH = F

#source("code/analysis/prepare.R")
#load("cases/panel-prepped.RData")
load("cases/panel-prepped_MLI.RData")

# select only lowest level
df <- df[df$lowest_level == 1, ]

# select only observations for month X
df$month <- months(strptime(df$Date, format='%Y-%m-%d'))
#df <- df[df$month == 'March', ]

# remove observations with the same weather
#df2 <- df %>% group_by(regid) %>% summarize(index=sum(t2m * tp)) %>% ungroup
#duplicates <- df2[duplicated(df2$index), 'regid']
#df <- df[!(df$regid %in% duplicates$regid), ]

## Models to report:
## instantaneous, no policy
mod.in <- felm(dlog ~ absh + de + t2m + tp + utci + utci2 + ssrd | factor(regid) + factor(regid) : factor(week) + factor(superset) : factor(Date) | 0 | regid, data=df)

#mod.in <- felm(dlog ~ absh + absh2 + de + de2 + q + q2 + r + r2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 | factor(regid) + factor(regid) : factor(week) + factor(superset) : factor(Date) | 0 | regid, data=df)
#mod.in <- felm(dlog ~ absh + de + q + r + ssrd + t2m + tp | factor(regid) + factor(regid) : factor(week) + factor(superset) : factor(Date) | 0 | regid, data=df)
summary(mod.in)

# identify problematic variable
#pairs(~absh + de + q + r + ssrd + t2m + tp, data=df)

## xval-delay, no policy
df.pred <- df
#weathers <- c('absh', 'absh2', 'de', 'de2', 'q', 'q2', 'r', 'r2', 'ssrd', 'ssrd2', 't2m', 't2m2', 'tp', 'tp2')
#for (weather in weathers)
#    df.pred[, weather] <- df[, paste0(weather, ".pred")]

if (SWITCH == T) {
weathers <- c('absh', 'r', 'ssrd', 't2m', 'tp', 'utci')
weathers_linearA <- c()
weathers_linearB <- c()
weathers_orders12A <- c()
weathers_orders12B <- c()
for (weather in weathers) {
	weathers_linearA <- c(weathers_linearA, paste0(weather, '.predA'))
	weathers_linearB <- c(weathers_linearB, paste0(weather, '.predB'))
	weathers_orders12A <- c(weathers_orders12A, paste0(weather, '.predA'), paste0(weather, '2.predA'))
	weathers_orders12B <- c(weathers_orders12B, paste0(weather, '.predB'), paste0(weather, '2.predB'))
}

formula <- 'dlog ~'
formulaA <- 'dlog ~'
formulaB <- 'dlog ~'
for (variable in weathers_orders12A) {
	formula <- paste0(formula, ' + ', variable)
	formulaA <- paste0(formulaA, ' + ', variable)
}
for (variable in weathers_orders12B) {
	formula <- paste0(formula, ' + ', variable)
	formulaB <- paste0(formulaB, ' + ', variable)
}
fe <- ' | factor(regid) + factor(regid) : factor(week) + factor(superset) : factor(Date) | 0 | regid'
formula <- paste0(formula, fe)
formulaA <- paste0(formulaA, fe)
formulaB <- paste0(formulaB, fe)

formula = as.formula(strwrap(formula, width=10000, simplify=T))
formulaA = as.formula(strwrap(formulaA, width=10000, simplify=T))
formulaB = as.formula(strwrap(formulaB, width=10000, simplify=T))

mod.xn <- felm(formula, data=df.pred)
summary(mod.xn)

mod.xnA <- felm(formulaA, data=df.pred)
summary(mod.xnA)

mod.xnB <- felm(formulaB, data=df.pred)
summary(mod.xnB)

library(stargazer)
stargazer(list(mod.xn, mod.xnA, mod.xnB), omit.stat='ser', no.space=T, column.labels=c("Both lags", "Lag 2-8", "Lag 10-20"))
}

## new on 2020-06-23

# select weather variables
weather_subset <- c('absh.predA', 'r.predA', 'ssrd.predA', 'tp.predA', 't2m.predA', 'utci.predA', 'utci2.predA',
					'absh.predB', 'r.predB', 'ssrd.predB', 'tp.predB', 't2m.predB', 'utci.predB', 'utci2.predB')

df[, c(weather_subset)] <- scale(df[, c(weather_subset)])

# read in policy fixed effects
poldf <- read.csv("policy/policyml.csv")
df <- merge(df, poldf, by.x=c('Country', 'Region', 'Locality', 'Date'), by.y=c('Country', 'Region', 'Locality', 'Date'))

# select observations
rows <- df$fex != '' & rowSums(is.na(df[, c('dlog', weather_subset)])) == 0

# project out fixed effects
dmdf <- demeanlist(df[rows, c('dlog', weather_subset)], list(factor(df$regid[rows]), factor(paste(df$regid[rows], df$week[rows])), factor(paste(df$superset[rows], df$Date[rows]))))

# estimate without policy FE
mod.AB <- lm(dlog ~ 0 + absh.predA + r.predA + ssrd.predA + tp.predA + t2m.predA + utci.predA + utci2.predA +
	absh.predB + r.predB + ssrd.predB + tp.predB + t2m.predB + utci.predB + utci2.predB, data=dmdf)
summary(mod.AB)

# estimate with policy FE
mod.ABpolicy <- felm(dlog ~ absh.predA + r.predA + ssrd.predA + tp.predA + t2m.predA + utci.predA + utci2.predA +
	absh.predB + r.predB + ssrd.predB + tp.predB + t2m.predB + utci.predB + utci2.predB | factor(regid) + factor(fex) + factor(superset) : factor(Date) | 0 | regid, data=df[rows, ])
summary(mod.ABpolicy)

tablat <- stargazer(list(mod.AB, mod.ABpolicy), omit.stat='ser', no.space=T, column.labels=c("standard FE", "policy FE"))
writeLines(tablat, "results/results1.tex")
