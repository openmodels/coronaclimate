library(dplyr)
library(lfe)
library(MASS)
library(glmnet)
library(stargazer)
library(boot)
library(sandwich)

## ===============================================================================
## no mobility, all observations
## ===============================================================================

DATAPATH <- './data/'
RESULTPATH <- './results/'
LABEL <- 'firstdiff_contemp_mobility'

## ===============================================================================

df <- read.csv(paste0(DATAPATH, "panel_prepared_firstdiff_contemp.csv"))

## only lowest level
df <- df[df$lowest_level == 1, ]
#df <- df[df$Country == 'USA', ]

## make day-of-week variable
df['dow'] <- df['days'] %% 7

# create vector of weather variables
weather.all <- c('absh', 'r', 't2m', 'ssrd', 'tp', 'utci')

# to get rid of bins
splitend <- function(x) strsplit(x, '_')[[1]][1]

# cost function for cross-validation
cost <- function(y, yhat=0) mean((y-yhat)**2.)

## prediction with user-specified cov-matrix
predict.rob <- function(x, clcov, newdata){
	if(missing(newdata)){ newdata <- x$model }
	tt <- terms(x)
	Terms <- delete.response(tt)
	m.mat <- model.matrix(Terms, data=newdata)
	m.coef <- x$coef
	fit <- as.vector(m.mat %*% x$coef)
	se.fit <- sqrt(diag(m.mat%*%clcov%*%t(m.mat)))
	return(list(fit=fit, se.fit=se.fit))
}

## START
df.variables <- read.csv('modelselection.csv')

mobvars <- c('mobility_pca1', 'mobility_pca2', 'mobility_pca3')

for (mobvar in mobvars) {
	print(mobvar)

	rows <- (rowSums(is.na(df[, c(mobvar, weather.all)])) == 0)
	vars <- c(mobvar, weather.all)

	# project out FE for fitting with the lm package further below (for predictions)
	if (nrow(df[rows, ]) < 100) {
		next
	}
	dmdf <- demeanlist(df[rows, vars], list(factor(paste(df$regid[rows], df$dow[rows])), factor(paste(df$regid[rows], df$week[rows])), factor(paste(df$superset[rows], df$Date[rows]))))

	for (specification in c('pol1')) {
		print(specification)
		models <- list()

		steps <- c(10)

		for (i.step in steps) {

			variables <- c()

			weather.step <- unlist(strsplit(as.character(df.variables[i.step,'variables']), ','))
			print(weather.step)

			if (specification == 'pol1') {
				variables <- c(variables, weather.step)
			}

			modelstring <- paste0(mobvar, paste('~ 0 +',  paste(variables, collapse=' + '), sep=''))

			# fit with felm package, and save coefficients, standard errors, and a few statistics in csv files
			modelformula <- as.formula(paste0(modelstring, ' | factor(regid) : factor(dow) + factor(regid) : factor(week) + factor(superset) : factor(Date) | 0 | regid'))
			print(modelformula)
			mod.felm <- felm(modelformula, data=df[rows, ], na.action=na.omit)
			s <- summary(mod.felm)
			write.csv(data.frame(s$coefficients), file=paste0(RESULTPATH, 'coeffs_', LABEL, '_', specification, '_', mobvar, '_', i.step, '.csv'))
			df.stats <- data.frame(index=1)
			df.stats["r2"] <- s$r2
			df.stats["r2adj"] <- s$r2adj
			df.stats["nobs"] <- s$N
			df.stats["ncoeffs"] <- s$p
			df.stats["ncoeffswithin"] <- s$Pp
			df.stats["SST"] <- sum((df[rows, mobvar] - mean(df[rows, mobvar]))^2.)
			df.stats["SSE"] <- sum(mod.felm$residuals^2.)
			df.stats["fstat"] <- s$fstat

			# get 10-fold crossvalidation error
			glm.fit <- glm(as.formula(modelstring), family='gaussian', data=dmdf)
			delta <- cv.glm(dmdf, glm.fit, cost, 10)$delta
			df.stats["delta1"] <- delta[1]
			df.stats["delta2"] <- delta[2]

			# fit with lm package, for r2-within and also for predictions below
			mod.lm <- lm(as.formula(modelstring), data=dmdf, na.action=na.omit)
			vcov.cluster <- vcovCL(mod.lm, dmdf$regid)
			df.stats["r2within"] <- summary(mod.lm)$r.squared
			df.stats["r2adjwithin"] <- summary(mod.lm)$adj.r.squared
			df.stats["SSTwithin"] <- sum((dmdf[, mobvar] - mean(dmdf[, mobvar]))^2.)

			write.csv(df.stats, file=paste0(RESULTPATH, 'stats_', LABEL, '_', specification, '_', mobvar, '_', i.step, '.csv'))
		}
	}
}