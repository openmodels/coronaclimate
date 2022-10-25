library(dplyr)
library(lfe)
#library(ggplot2)
library(MASS)
library(glmnet)
library(stargazer)
library(boot)
library(sandwich)

## ===============================================================================
## no mobility, all observations
## ===============================================================================

LAGS <- 1:10
DATAPATH <- './data/'
RESULTPATH <- './results/'
LABEL <- 'firstdiff_mortality_L10'

## ===============================================================================

df <- read.csv(paste0(DATAPATH, "panel_prepared_distlags_firstdiff_L10.csv"))

## only lowest level
df <- df[df$lowest_level == 1, ]

## make day-of-week variable
df['dow'] <- df['days'] %% 7

# create vector of weather variables
weather.all <- c('absh', 'r', 't2m', 'ssrd', 'tp', 'utci')

weather.lags <- c()
for (weather in weather.all) {
	for (ilag in LAGS) {
		weather.lags <- c(weather.lags, paste0(weather, '_L', ilag))
	}
}

weather.lags.sq <- c()
for (weather in weather.all) {
	for (ilag in LAGS) {
		weather.lags.sq <- c(weather.lags.sq, paste0(weather, '2_L', ilag))
		df[paste0(weather, '2_L', ilag)] <- df[paste0(weather, '_L', ilag)]**2.
	}
}


weather_lags <- function(weather.vector, lags) {
	weather.lags <- c()
	for (weather in weather.vector) {
		for (ilag in lags) {
			weather.lags <- c(weather.lags, paste0(weather, '_L', ilag))
		}
	}
	weather.lags
}

weather_lags_sq <- function(weather.vector, lags) {
	weather.lags <- c()
	for (weather in weather.vector) {
		for (ilag in lags) {
			weather.lags <- c(weather.lags, paste0(weather, '2_L', ilag))
		}
	}
	weather.lags
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
    	coeff1 <- coeff1 + coef[paste0(weather, '_L', lag), 1]
    	var1 <- var1 + vcovM[paste0(weather, '_L', lag), paste0(weather, '_L', lag)]
	    for (other.lag in (1:nlags)[(1:nlags) > lag]) {
	      if (other.lag >= 1) {
	      	var1 <- var1 + 2. * vcovM[paste0(weather, '_L', lag), paste0(weather, '_L', other.lag)]
	      }
	    }
	}
	se1 <- sqrt(var1)
	return(list(coeffs=c(coeff1), se=c(se1)))
}

cumulative_beta_pol2 <- function(modelobj, vcovM, weather, nlags){
	coef <- summary(modelobj)$coefficients
	coeff1 <- 0.
	coeff2 <- 0.
	var1 <- 0.
	var2 <- 0.
	for (lag in 1:nlags) {
    	coeff1 <- coeff1 + coef[paste0(weather, '_L', lag), 1]
    	coeff2 <- coeff2 + coef[paste0(weather, '2_L', lag), 1]
    	var1 <- var1 + vcovM[paste0(weather, '_L', lag), paste0(weather, '_L', lag)]
    	var2 <- var2 + vcovM[paste0(weather, '2_L', lag), paste0(weather, '2_L', lag)]
	    for (other.lag in (1:nlags)[(1:nlags) > lag]) {
	      if (other.lag >= 1) {
	      	var1 <- var1 + 2. * vcovM[paste0(weather, '_L', lag), paste0(weather, '_L', other.lag)]
	      	var2 <- var2 + 2. * vcovM[paste0(weather, '2_L', lag), paste0(weather, '2_L', other.lag)]
	      }
	    }
	}
	se1 <- sqrt(var1)
	se2 <- sqrt(var2)
	return(list(coeffs=c(coeff1, coeff2), se=c(se1, se2)))
}

## demean and scale, only for polynomical models
df.norm <- data.frame(matrix(ncol=length(weather.lags), nrow=2))
colnames(df.norm) <- weather.lags
rownames(df.norm) <- c('mean', 'std')

for (weather in weather.lags) {
	df.norm[1, weather] <- mean(df[, weather], na.rm=TRUE)
	df.norm[2, weather] <- sd(df[, weather], na.rm=TRUE)
}

write.csv(df.norm, paste0(RESULTPATH, 'means-sds_', LABEL, '.csv'))

df[, c(weather.lags)] <- scale(df[, c(weather.lags)])

## START
df.variables <- read.csv('modelselection.csv')

mobmodes <- c('nomob')

for (mobmode in mobmodes) {
	print(mobmode)
	if (mobmode == 'nomob') {
		vars <- c('mort.4w', weather.lags, weather.lags.sq)
		rows <- rowSums(is.na(df[, vars])) == 0

	} else if (mobmode == 'withoutmob') {
		rows <- (rowSums(is.na(df[, c('mort.4w', weather.lags)])) == 0) & (rowSums(is.na(df[, mobility.lags])) != 0)
		vars <- c('mort.4w', weather.lags, weather.lags.sq, mobility.lags)

	} else if (mobmode %in% c('mobcntr', 'moball', 'mobpca')) {
		rows <- rowSums(is.na(df[, c('mort.4w', weather.lags, mobility.lags)])) == 0
		vars <- c('mort.4w', weather.lags, weather.lags.sq, mobility.lags)
	}

	# calculate means of all variables and store in dataframe
	df.means <- data.frame(matrix(ncol=length(c(weather.lags, weather.lags.sq)), nrow=1))
	colnames(df.means) <- c(weather.lags, weather.lags.sq)
	for (weather in c(weather.lags, weather.lags.sq)) {
		df.means[1, weather] <- mean(df[rows, weather], na.rm=TRUE)
	}

	# project out FE for fitting with the lm package further below (for predictions)
	if (nrow(df[rows, ]) < 100) {
		next
	}
	dmdf <- demeanlist(df[rows, vars], list(factor(paste(df$regid[rows], df$dow[rows])), factor(paste(df$regid[rows], df$week[rows])), factor(paste(df$superset[rows], df$Date[rows]))))

	for (specification in c('pol1')) { #, 'pol2'
		print(specification)
		models <- list()

		steps <- c(10)

		for (i.step in steps) {
			variables <- c()

			weather.step <- unlist(strsplit(as.character(df.variables[i.step,'variables']), ','))
			print(weather.step)

			if (mobmode %in% c('nomob', 'withoutmob', 'mobcntr')) {
				variables <- c(variables)
			} else if (mobmode == 'moball') {
				variables <- c(variables, mobility.cat.lags)
			} else if (mobmode == 'mobpca') {
				variables <- c(variables, mobility.pca.lags)
			}

			if (specification == 'pol1') {
				variables <- c(variables, weather_lags(weather.step, LAGS))
			} else if (specification == 'pol2') {
				variables <- c(variables, weather_lags(weather.step, LAGS), weather_lags_sq(weather.step, LAGS))
			}

			modelstring <- paste('mort.4w ~ 0 +',  paste(variables, collapse=' + '), sep='')

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
			df.stats["ncoeffs"] <- s$p
			df.stats["ncoeffswithin"] <- s$Pp
			df.stats["SST"] <- sum((df[rows, 'mort.4w'] - mean(df[rows, 'mort.4w']))^2.)
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
			df.stats["SSTwithin"] <- sum((dmdf[, 'mort.4w'] - mean(dmdf[, 'mort.4w']))^2.)
			write.csv(df.stats, file=paste0(RESULTPATH, 'stats_', LABEL, '_', specification, '_', mobmode, '_', i.step, '.csv'))

			# save predicted values in csv files
			for (lag in c(LAGS, 'cum5', 'cum')) {

				for (weather in weather.step) {

					weather.lag <- paste0(weather, '_L1')
					weather.lag.sq <- paste0(weather, '2_L1')

					new.data <- data.frame(matrix(ncol=length(variables), nrow=100))
					colnames(new.data) <- variables
					new.data[,] <- 0.

					coords_scaled <- seq(from=min(df[, weather.lag], na.rm=TRUE), to=max(df[, weather.lag], na.rm=TRUE), length=100)
					if (specification == 'pol1') {
						coords_scaled_demeaned <- coords_scaled - df.means[1, weather.lag]
						coords_unscaled <- coords_scaled * df.norm[2, weather.lag] + df.norm[1, weather.lag]
					} else if (specification == 'pol2') {
						coords_scaled_demeaned <- coords_scaled - df.means[1, weather.lag]
						coords_scaled_demeaned_sq <- (coords_scaled**2.) - df.means[1, weather.lag.sq]
						coords_unscaled <- coords_scaled * df.norm[2, weather.lag] + df.norm[1, weather.lag]
					}

					if (lag == 'cum') { 
						for (ilag in LAGS) {
							if (specification == 'pol1') {
								new.data[, paste0(weather, '_L', ilag)] <- coords_scaled_demeaned
							} else if (specification == 'pol2') {
								new.data[, paste0(weather, '_L', ilag)] <- coords_scaled_demeaned
								new.data[, paste0(weather, '2_L', ilag)] <- coords_scaled_demeaned_sq	
							}
						}
						if (specification == 'pol1') {
							retvals <- cumulative_beta_pol1(mod.lm, vcov.cluster, weather, length(LAGS))
						} else {
							retvals <- cumulative_beta_pol2(mod.lm, vcov.cluster, weather, length(LAGS))
						}
						df.cum <- data.frame(retvals)
						write.csv(df.cum, paste0(RESULTPATH, 'coeffscum_', LABEL, '_', specification, '_L-', lag, '_', mobmode, '_', i.step, '_', weather, '.csv'))
					} else if (lag == 'cum5') { 
						for (ilag in 1:5) {
							if (specification == 'pol1') {
								new.data[, paste0(weather, '_L', ilag)] <- coords_scaled_demeaned
							} else if (specification == 'pol2') {
								new.data[, paste0(weather, '_L', ilag)] <- coords_scaled_demeaned
								new.data[, paste0(weather, '2_L', ilag)] <- coords_scaled_demeaned_sq	
							}
						}
						if (specification == 'pol1') {
							retvals <- cumulative_beta_pol1(mod.lm, vcov.cluster, weather, 5)
						} else {
							retvals <- cumulative_beta_pol2(mod.lm, vcov.cluster, weather, 5)
						}
						df.cum <- data.frame(retvals)
						write.csv(df.cum, paste0(RESULTPATH, 'coeffscum_', LABEL, '_', specification, '_L-', lag, '_', mobmode, '_', i.step, '_', weather, '.csv'))
					} else {
						if (specification == 'pol1') {
							new.data[, paste0(weather, '_L', lag)] <- coords_scaled_demeaned
						} else if (specification == 'pol2') {
							new.data[, paste0(weather, '_L', lag)] <- coords_scaled_demeaned
							new.data[, paste0(weather, '2_L', lag)] <- coords_scaled_demeaned_sq	
						}
					}
					df.predict <- data.frame(predict.rob(mod.lm, vcov.cluster, new.data))
					df.predict[, weather.lag] <- coords_unscaled
					write.csv(df.predict, paste0(RESULTPATH, 'predict_', LABEL, '_', specification, '_L-', lag, '_', mobmode, '_', i.step, '_', weather, '.csv'))
				}
			}

			# save the residuals as csv
			df.residuals <- data.frame(residuals=residuals.lm(mod.lm, type='response'))
			df.residuals$GEO_ID <- df[rows, "GEO_ID"]
			df.residuals$GEO_ID_REFERENCE <- df[rows, "GEO_ID_REFERENCE"]
			df.residuals$days <- df[rows, "days"]
			write.csv(df.residuals, file=paste0(RESULTPATH, 'residuals_', LABEL, '_', specification, '_', mobmode, '_', i.step, '.csv'))
			models[[length(models)+1]] <- mod.felm

		}
	}
}