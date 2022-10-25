library(dplyr)
library(lfe)
#library(ggplot2)
library(MASS)
library(glmnet)
library(stargazer)
library(boot)
library(sandwich) # for robust standard errors

## ===============================================================================
## no mobility, all observations
## ===============================================================================

LAGS <- 1:10
NBINS <- 5
DATAPATH <- './data/'
RESULTPATH <- './results/'
LABEL <- 'distlag_bins_L10'

## ===============================================================================

df <- read.csv(paste0(DATAPATH, "panel_prepared_distlags_bins_firstdiff_L10.csv"))

## only lowest level
df <- df[df$lowest_level == 1, ]
#df <- df[df$Country == 'USA', ]

## make day-of-week variable
df['dow'] <- df['days'] %% 7

# create vector of weather variables
weather.all <- c('absh', 'r', 't2m', 'ssrd', 'tp', 'utci')

print('making binned variables')
weather.variables <- c()
for (weather in weather.all) {
	for (ilag in LAGS) {
		nbins <- NBINS
		boundaries <- seq(from=0, to=1, length=nbins+1)
		for (i.bin in 1:nbins) {
			weather.variables <- c(weather.variables, paste0(weather, '_', i.bin, '_L', ilag))	
		}
	}
}


weather_binned_lags <- function(weather.vector, lags, nbins) {
	weather.variables <- c()
	for (weather in weather.vector) {
		for (ilag in lags) {
			boundaries <- seq(from=0, to=1, length=nbins+1)
			for (i.bin in 1:nbins) {
				if (i.bin != 3) {
					weather.variables <- c(weather.variables, paste0(weather, '_', i.bin, '_L', ilag))
				}
			}
		}
	}
	weather.variables
}

weather_bins <- function(weather.vector, nbins) {
	weather.variables <- c()
	for (weather in weather.vector) {
		for (i.bin in 1:nbins) {
			if (i.bin != 3) {
				weather.variables <- c(weather.variables, paste0(weather, '_', i.bin))
			}
		}
	}
	weather.variables
}


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

cumulative_beta_pol1 <- function(modelobj, vcovM, weather, nlags){
	coef <- summary(modelobj)$coefficients
	coeff1 <- 0.
	var1 <- 0.
	for (lag in 1:nlags) {
		#print(paste0(weather, '_L', lag))
		if (!(paste0(weather, '_L', lag) %in% labels(coef)[[1]])) {
			coeff1 <- NA
			break
		}
    	coeff1 <- coeff1 + coef[paste0(weather, '_L', lag), 1]
    	var1 <- var1 + vcovM[paste0(weather, '_L', lag), paste0(weather, '_L', lag)]
	    for (other.lag in (1:nlags)[(1:nlags) > lag]) {
	      if (other.lag >= 1) {
		     if (!(paste0(weather, '_L', other.lag) %in% labels(vcovM)[[2]])) {
				var1 <- NA
				break
			}
	      	var1 <- var1 + 2. * vcovM[paste0(weather, '_L', lag), paste0(weather, '_L', other.lag)]
	      }
	    }
	}
	se1 <- sqrt(var1)
	return(list(coeffs=c(coeff1), se=c(se1)))
}

mobility.pca.lags <- c()
for (mobility_var in c('mobility_pca1')) {
	for (ilag in LAGS) {
		mobility.pca.lags <- c(mobility.pca.lags, paste0(mobility_var, '_L', ilag))
	}
}

mobility.cat.lags <- c()
for (mobility_var in c('mobility_retail_and_recreation', 'mobility_grocery_and_pharmacy', 'mobility_parks',
				'mobility_transit_stations', 'mobility_workplaces', 'mobility_residential')) {
	for (ilag in LAGS) {
		mobility.cat.lags <- c(mobility.cat.lags, paste0(mobility_var, '_L', ilag))
	}
}

mobility.lags <- c(mobility.pca.lags, mobility.cat.lags)

dfc <- df # save a data frame that remains constant for use in several loops below

# START
df.variables <- read.csv('modelselection.csv')

mobmodes <- c('nomob', 'withoutmob', 'mobcntr', 'moball', 'mobpca')

for (mobmode in mobmodes) {
	df <- dfc
	
	print(mobmode)
	if (mobmode == 'nomob') {
		vars <- c('dlog', weather.variables)
		rows <- rowSums(is.na(df[, vars])) == 0

	} else if (mobmode == 'withoutmob') {
		rows <- (rowSums(is.na(df[, c('dlog', weather.variables)])) == 0) & (rowSums(is.na(df[, mobility.lags])) != 0)
		vars <- c('dlog', weather.variables, mobility.lags)

	} else if (mobmode %in% c('mobcntr', 'moball', 'mobpca')) {
		rows <- rowSums(is.na(df[, c('dlog', weather.variables, mobility.lags)])) == 0
		vars <- c('dlog', weather.variables, mobility.lags)
	}
	# project out FE for fitting with the lm package further below (for predictions)
	if (nrow(df[rows, ]) < 100) {
		next
	}
	dmdf <- demeanlist(df[rows, vars], list(factor(paste(df$regid[rows], df$dow[rows])), factor(paste(df$regid[rows], df$week[rows])), factor(paste(df$superset[rows], df$Date[rows]))))

	for (specification in c('bins')) {
		print(specification)
		models <- list()

		steps <- unique(df.variables$n)
		steps <- c(10)

		for (i.step in steps) {
			variables <- c()

			weather.step <- unlist(strsplit(as.character(df.variables[df.variables$n == i.step, 'variables']), ','))
			print(weather.step)

			if (mobmode %in% c('nomob', 'withoutmob', 'mobcntr')) {
				variables <- c(variables)
			} else if (mobmode == 'moball') {
				variables <- c(variables, mobility.cat.lags)
			} else if (mobmode == 'mobpca') {
				variables <- c(variables, mobility.pca.lags)
			}

			if (specification == 'bins') {
				variables <- c(variables, weather_binned_lags(weather.step, LAGS, NBINS))
			}

			modelstring <- paste('dlog ~ 0 +',  paste(variables, collapse=' + '), sep='')

			# fit with felm package, and save coefficients, standard errors, and a few statistics in csv files
			modelformula <- as.formula(paste0(modelstring, ' | factor(regid) : factor(dow) + factor(regid) : factor(week) + factor(superset) : factor(Date) | 0 | regid'))
			print(modelformula)
			mod.felm <- felm(modelformula, data=df[rows, ], na.action=na.omit)
			s <- summary(mod.felm)
			write.csv(data.frame(s$coefficients), file=paste0(RESULTPATH, 'coeffs_', LABEL, '_', specification, '_', mobmode, '_', i.step, '.csv'))
			df.stats <- data.frame(index=1)
			df.stats["r2"] <- s$r2
			df.stats["r2adj"] <- s$r2adj
			df.stats["nobs"] <- s$N

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
			write.csv(df.stats, file=paste0(RESULTPATH, 'stats_', LABEL, '_', specification, '_', mobmode, '_', i.step, '.csv'))

			# save predicted values in csv files
			for (lag in c('cum5', 'cum')) {

				df.cum <- data.frame(matrix(nrow=NBINS, ncol=3))
				names(df.cum) <- c('variable', 'coef', 'se')
				for (weather in weather_bins(weather.step, NBINS)) {
					#print(weather)
					if (lag == 'cum') { 
						retvals <- cumulative_beta_pol1(mod.lm, vcov.cluster, weather, length(LAGS))
						df.cum <- rbind(df.cum, data.frame(variable=weather, coef=retvals$coeffs, se=retvals$se))

					} else if (lag == 'cum5') { 
						retvals <- cumulative_beta_pol1(mod.lm, vcov.cluster, weather, 5)
						df.cum <- rbind(df.cum, data.frame(variable=weather, coef=retvals$coeffs, se=retvals$se))
					}
				}
				write.csv(df.cum, paste0(RESULTPATH, 'coeffscum_', LABEL, '_', specification, '_L-', lag, '_', mobmode, '_', i.step, '.csv'))
			}
		}
	}
}