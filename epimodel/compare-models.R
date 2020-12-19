setwd("~/Dropbox/Coronavirus and Climate/code/epimodel")

library(lubridate)

result.files <- c('epimodel-0604.csv', 'epimodel-0604-ww.csv', 'epimodel-0615.csv',
                  'epimodel-meta-0616.csv', 'epimodel-0803.csv', 'epimodel-meta-0806-x5.csv',
                  'epimodel-meta-0817.csv', 'epimodel-meta-0817-noprior.csv',
                  'epimodel-meta-0907-nobs.csv', 'epimodel-meta-0907-pop.csv', 'epimodel-meta-0907-region.csv',
                  'epimodel-meta-0921-pop.csv', 'epimodel-meta-1030-all-pop.csv', 'epimodel-meta-1111-mixed-all-pop.csv')
result.names <- c('Kucharski et al.', 'Include Weather', 'Handle Deaths',
                  'Death & Not', 'Early/Late & Mobility', 'Const. Omega & Priors',
                  'New Weather', 'Drop Priors',
                  'Omega Effect, by Pop.', 'Omega Effect, by Obs.', 'Omega Effect, by Reg.',
                  'Smooth Omega', 'Estimated Testing', 'OLS Compare')

coeff.names <- list('alpha'='Time Variance', 'invsigma'='Incubation Period', 'invgamma'='Infectious Period', 'omega'='Reporting Rate', 'error'='Model Error', 'e.absh'='Beta[Abs. Humid.]', 'e.r'='Beta[Rel. Humid.]', 'e.t2m'='Beta[Surface Temp.]', 'e.tp'='Beta[Total Prec.]', 'o.absh'='Zeta[Abs. Humid.]', 'o.r'='Zeta[Rel. Humid.]', 'o.t2m'='Zeta[Surface Temp.]', 'o.tp'='Zeta[Total Prec.]', 'deathrate'='Death Rate', 'deathomegaplus'='Death Reporting', 'portion_early'='Portion Early', 'mobility_slope'='Mobility Slope', 'logbeta'='Log Transmission', 'logomega'='Log Reporting', 'eein'='Infected Imports', 'omega0'='Initial Reporting', 'domega'='Reporting Slope', 'invkappa'='Testing Delay', 'invtheta'='Reporting Delay', 'e.ssrd', 'e.ssrd'='Beta[Solar Rad.]', 'e.utci'='Beta[Thermal Index]', 'o.ssrd'='Zeta[Solar Rad.]', 'o.utci'='Zeta[Thermal Index]')
coeff.order <- c('Time Variance', 'Incubation Period', 'Infectious Period', 'Testing Delay', 'Reporting Delay', 'Reporting Rate', 'Initial Reporting', 'Reporting Slope', 'Model Error', 'Beta[Abs. Humid.]', 'Beta[Rel. Humid.]', 'Beta[Surface Temp.]', 'Beta[Total Prec.]', 'Beta[Solar Rad.]', 'Beta[Thermal Index]', 'Zeta[Abs. Humid.]', 'Zeta[Rel. Humid.]', 'Zeta[Surface Temp.]', 'Zeta[Total Prec.]', 'Zeta[Solar Rad.]', 'Zeta[Thermal Index]', 'Death Rate', 'Death Reporting', 'Portion Early', 'Mobility Slope', 'Log Transmission', 'Log Reporting', 'Infected Imports')

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

        if (rr >= which(result.files == 'epimodel-0803.csv') & rr <= which(result.files == 'epimodel-meta-0806-x5.csv')) {
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
    performs <- rbind(performs, data.frame(model=result.names[rr], regions, rhat))
}

performs

library(ggplot2)

fits$model <- factor(fits$model, levels=rev(result.names))
fits$param.name <- factor(fits$param.name, coeff.order)

gp <- ggplot(fits[!(fits$param %in% c('logbeta', 'logomega')),], aes(model, meanmu)) +
    facet_wrap(~ param.name, scales='free_x', ncol=8) +
    geom_errorbar(aes(ymin=mean25, ymax=mean75)) +
    geom_linerange(aes(ymin=mu25q, ymax=mu75q), col='red', linetype='dashed') +
    geom_linerange(aes(ymin=ymin, ymax=ymax), col='blue') +
    coord_flip() + scale_y_continuous(expand=c(.01, .01)) + theme_bw() +
    xlab(NULL) + ylab("Parameter Estimate")
ggsave("../../figures/compare-models.pdf", gp, width=10, height=8)
