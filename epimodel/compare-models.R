setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

library(lubridate)

result.files <- c('epimodel-0604.csv', 'epimodel-0604-ww.csv', 'epimodel-0614.csv', 'epimodel-0615.csv', 'epimodel-0616.csv', 'epimodel-meta-0616.csv', 'epimodel-0728.csv', 'epimodel-0730.csv', 'epimodel-0803.csv', 'epimodel-0805.csv', 'epimodel-0806.csv', 'epimodel-0806-x5.csv', 'epimodel-meta-0806-x5.csv', 'epimodel-0817-run1.csv', 'epimodel-meta-0817-run1.csv', 'epimodel-0817.csv', 'epimodel-meta-0817.csv', 'epimodel-0817-noprior.csv', 'epimodel-meta-0817-noprior.csv', 'epimodel-0907.csv', 'epimodel-meta-0907.csv')
result.names <- c('Kucharski et al.', 'Include Weather', 'Handle Deaths', 'Improve Priors', 'Death & Not', 'Full Meta', 'Early & Late', 'Drop Reportable', 'Mobility Prior', 'Constant Omega', 'Weather Priors', 'Loose Priors', 'Full Meta 2', 'New Weather', 'Full Meta 3', 'New Weather R2', 'Full Meta 3 R2', 'New Weather NoP', 'Full Meta 3 NoP', 'Omega Effect', 'Full Meta 4')
result.ctime <- parse_date_time(c('06/04/2020 20:05:51', '06/11/2020 07:49:16', '06/14/2020 14:03:17', '06/15/2020 17:26:33', '06/17/2020 17:59:03', '06/18/2020 15:17:53', '07/30/2020 13:09:45', '08/01/2020 10:52:20', '08/03/2020 18:33:48', '08/06/2020 00:05:45', '08/06/2020 11:21:01', '08/06/2020 23:03:18', '08/13/2020 13:16:36', NA, NA, NA, NA, NA, NA, NA, NA), '%m/%d/%y %H:%M:%S', tz='GB')

coeff.names <- list('alpha'='Time Variance', 'invsigma'='Incubation Period', 'invgamma'='Infectious Period', 'invkappa'='Reporting Delay', 'omega'='Reporting Rate', 'error'='Model Error', 'e.absh'='Beta[Abs. Humid.]', 'e.r'='Beta[Rel. Humid.]', 'e.t2m'='Beta[Surface Temp.]', 'e.tp'='Beta[Total Prec.]', 'o.absh'='Zeta[Abs. Humid.]', 'o.r'='Zeta[Rel. Humid.]', 'o.t2m'='Zeta[Surface Temp.]', 'o.tp'='Zeta[Total Prec.]', 'deathrate'='Death Rate', 'deathomegaplus'='Death Reporting', 'portion_early'='Portion Early', 'mobility_slope'='Mobility Slope', 'logbeta'='Log Transmission', 'eein'='Infected Imports')
coeff.order <- c('Time Variance', 'Incubation Period', 'Infectious Period', 'Reporting Delay', 'Reporting Rate', 'Model Error', 'Beta[Abs. Humid.]', 'Beta[Rel. Humid.]', 'Beta[Surface Temp.]', 'Beta[Total Prec.]', 'Zeta[Abs. Humid.]', 'Zeta[Rel. Humid.]', 'Zeta[Surface Temp.]', 'Zeta[Total Prec.]', 'Death Rate', 'Death Reporting', 'Portion Early', 'Mobility Slope', 'Log Transmission', 'Infected Imports')

fits <- data.frame()
performs <- data.frame()
for (rr in 1:length(result.files)) {
    filepath <- file.path("../../results", result.files[rr])
    print(filepath)
    df <- read.csv(filepath)
    for (param in unique(df$param)) {
        if ("  " %in% df$regid[df$param == param]) {
            meanmu <- df$mu[df$param == param & df$regid == "  "]
            mean25 <- df$ci25[df$param == param & df$regid == "  "]
            mean75 <- df$ci75[df$param == param & df$regid == "  "]
        } else {
            meanmu <- mean(df$mu[df$param == param], na.rm=T)
            mean25 <- mean(df$ci25[df$param == param], na.rm=T)
            mean75 <- mean(df$ci75[df$param == param], na.rm=T)
        }
        mu25q <- quantile(df$mu[df$param == param], .25, na.rm=T)
        mu75q <- quantile(df$mu[df$param == param], .75, na.rm=T)
        if ('rhat' %in% names(df))
            rhat <- mean(df$rhat[df$param == param], na.rm=T)
        else
            rhat <- NA

        if (rr >= which(result.files == 'epimodel-0730.csv') & rr <= which(result.files == 'epimodel-0806-x5.csv')) {
            if (param == 'e.absh')
                param <- 'e.r'
            else if (param == 'e.r')
                param <- 'e.t2m'
        }

        ymin <- NA
        ymax <- NA
        if (substring(as.character(param), 1, 2) == 'e.')
            if (!is.na(mean25) && (mean75 > 1 || mean25 < -1 || mu75q > 1 || mu25q < -1)) {
                ymin <- -1
                ymax <- 1
                meanmu <- NA
                mean25 <- NA
                mean75 <- NA
                mu25q <- NA
                mu75q <- NA
            }

        fits <- rbind(fits, data.frame(model=result.names[rr], param, param.name=coeff.names[[param]], meanmu, mean25, mean75, mu25q, mu75q, ymin, ymax, rhat=rhat))
    }
    if ('rhat' %in% names(df))
        rhat <- median(df$rhat, na.rm=T)
    else
        rhat <- NA
    regions <- length(unique(df$regid))
    time0 <- result.ctime[rr]
    time1 <- file.info(filepath)$mtime
    rate <- (time1 - time0) / (regions - 1)
    performs <- rbind(performs, data.frame(model=result.names[rr], regions, rate, rhat))
}

performs[, -3]

library(ggplot2)

fits$model <- factor(fits$model, levels=rev(result.names))
fits$param.name <- factor(fits$param.name, coeff.order)

ggplot(fits, aes(model, meanmu)) +
    facet_wrap(~ param.name, scales='free_x', ncol=7) +
    geom_errorbar(aes(ymin=mean25, ymax=mean75)) +
    geom_linerange(aes(ymin=mu25q, ymax=mu75q), col='red', linetype='dashed') +
    geom_linerange(aes(ymin=ymin, ymax=ymax), col='blue') +
    coord_flip() + scale_y_continuous(expand=c(.01, .01)) + theme_bw() +
    xlab(NULL) + ylab("Parameter Estimate")

