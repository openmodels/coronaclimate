setwd("~/Dropbox/Coronavirus and Climate")

##source("code/analysis/prepare.R")
#load("cases/panel-prepped.RData")
load("cases/panel-prepped_MLI.RData")

df <- df[df$lowest_level == 1, ]
df <- df %>%
    mutate(month = months(strptime(Date, format='%Y-%m-%d')))

df <- df[df$month == 'March', ]

## Models to report:
## instantaneous, no policy
mod.in <- felm(dlog ~ absh + absh2 + de + de2 + q + q2 + r + r2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 | factor(regid) + factor(regid) : factor(week) + factor(superset) : factor(Date) | 0 | regid, data=df)
summary(mod.in)

## xval-delay, no policy
df.pred <- df
#weathers <- c('absh', 'absh2', 'de', 'de2', 'q', 'q2', 'r', 'r2', 'ssrd', 'ssrd2', 't2m', 't2m2', 'tp', 'tp2')
#for (weather in weathers)
#    df.pred[, weather] <- df[, paste0(weather, ".pred")]

weathers <- c('absh', 'de', 'q', 'r', 'ssrd', 't2m', 'tp')
weathers_linearA <- c()
weathers_linearB <- c()
for (weather in weathers) {
	weathers_linearA <- c(weathers_linearA, paste0(weather, '.predA'))
	weathers_linearB <- c(weathers_linearB, paste0(weather, '.predB'))
}

formula <- 'dlog ~'
formulaA <- 'dlog ~'
formulaB <- 'dlog ~'
for (variable in weathers_linearA) {
	formula <- paste0(formula, ' + ', variable)
	formulaA <- paste0(formulaA, ' + ', variable)
}
for (variable in weathers_linearB) {
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


## instantaneous, with policy
poldf <- read.csv("policy/policyml.csv")
df$fex <- poldf$fex
mod.ip <- felm(dlog ~ absh + absh2 + de + de2 + q + q2 + r + r2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 | factor(regid) + factor(fex) + factor(superset) : factor(Date) | 0 | regid, data=subset(df, fex != ''))
summary(mod.ip)

## xval-delay, with policy
df.pred$fex <- df$fex
mod.xp <- felm(dlog ~ absh + absh2 + de + de2 + q + q2 + r + r2 + ssrd + ssrd2 + t2m + t2m2 + tp + tp2 | factor(regid) + factor(fex) + factor(superset) : factor(Date) | 0 | regid, data=subset(df.pred, fex != ''))
summary(mod.xp)

library(stargazer)
stargazer(list(mod.in, mod.xn, mod.ip, mod.xp), omit.stat='ser', no.space=T, column.labels=c("Inst.", "X-Val", "Inst. P.", "X-Val P."), add.lines=list(c("F-statistic(full model)", 2.56, 2.547, 42.9, 39.17), c("F-statistic(proj model)", 0.9522, 1.611, 2.436, 2.145)))

## Look for best performing predictors

setwd("~/Dropbox/Coronavirus and Climate")

library(lfe)

##source("code/analysis/prepare.R")
load("cases/panel-prepped.RData")

poldf <- read.csv("policy/policyml.csv")
df$fex <- poldf$fex

dmdf1 <- df[, c('dlog', 'absh.pred', 'de.pred', 'q.pred',  'r.pred', 'ssrd.pred', 't2m.pred', 'tp.pred')]
rows <- df$fex != '' & rowSums(is.na(dmdf1)) == 0
dmdf2 <- demeanlist(dmdf1[rows,], list(factor(df$regid[rows]), factor(df$fex[rows]), factor(paste(df$superset[rows], df$Date[rows]))))

mod.lfe <- felm(dlog ~ absh.pred + de.pred + q.pred + r.pred + ssrd.pred + tp.pred + t2m.pred | factor(regid) + factor(fex) + factor(superset) : factor(Date) | 0 | regid, data=df[rows,])
summary(mod.lfe)

mod.all <- lm(dlog ~ 0 + absh.pred + de.pred + q.pred + r.pred + ssrd.pred + tp.pred + t2m.pred, data=dmdf2)

library(MASS)
stepAIC(mod2, direction='both')

mod.aic <- lm(dlog ~ 0 + absh.pred + q.pred + r.pred + ssrd.pred + tp.pred, data=dmdf2)

library(glmnet)

cv.out <- cv.glmnet(as.matrix(dmdf2[, -1]), dmdf2[, 1], alpha=1, nfolds=3)
plot(cv.out$nzero, cv.out$cvm)
plot(cv.out)

bestlam <- cv.out$lambda.min

lasso.mod <- glmnet(as.matrix(dmdf2[, -1]), dmdf2[, 1], alpha=1)

lasso.pred <- predict(lasso.mod, s=bestlam, newx=as.matrix(dmdf2[, -1]))
mean((lasso.pred-dmdf2[, 1])^2)

predict(lasso.mod, type='coefficients', s=bestlam)

mod.lso <- lm(dlog ~ 0 + absh.pred + de.pred + r.pred + ssrd.pred + t2m.pred + tp.pred, data=dmdf2)

library(stargazer)
stargazer(list(mod.all, mod.aic, mod.lso), omit.stat='ser', no.space=T, column.labels=c("All", "AIC", "LASSO"))

anv.all <- data.frame(coeff=row.names(anova(mod.all)), varexp=1e6 * anova(mod.all)[,2] / sum(anova(mod.all)[,2]))
anv.aic <- data.frame(coeff=row.names(anova(mod.aic)), varexp=1e6 * anova(mod.aic)[,2] / sum(anova(mod.aic)[,2]))
anv.lso <- data.frame(coeff=row.names(anova(mod.lso)), varexp=1e6 * anova(mod.lso)[,2] / sum(anova(mod.lso)[,2]))

library(dplyr)
anv <- anv.all %>% left_join(anv.aic, by='coeff', suffix=c('.all', '.aic')) %>%
    left_join(anv.lso, by='coeff', suffix=c('', '.lso'))
library(xtable)
print(xtable(anv[-nrow(anv),]), include.rownames=F)
