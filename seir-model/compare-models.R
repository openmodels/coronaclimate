setwd("~/research/coronavirus/code/seir-model")

library(lubridate)

## result.files <- c('epimodel-0604.csv', 'epimodel-0604-ww.csv', 'epimodel-0615.csv',
##                   'epimodel-meta-0616.csv', 'epimodel-0803.csv', 'epimodel-meta-0806-x5.csv',
##                   'epimodel-meta-0817.csv', 'epimodel-meta-0817-noprior.csv',
##                   'epimodel-meta-0907-nobs.csv', 'epimodel-meta-0907-pop.csv', 'epimodel-meta-0907-region.csv',
##                   'epimodel-meta-0921-pop.csv', 'epimodel-meta-1030-all-pop.csv', 'epimodel-meta-1111-mixed-all-pop.csv')
## result.names <- c('Kucharski et al.', 'Include Weather', 'Handle Deaths',
##                   'Death & Not', 'Early/Late & Mobility', 'Const. Omega & Priors',
##                   'New Weather', 'Drop Priors',
##                   'Omega Effect, by Pop.', 'Omega Effect, by Obs.', 'Omega Effect, by Reg.',
##                   'Smooth Omega', 'Estimated Testing', 'OLS Compare')
result.files <- c('epimodel-meta-0314-noweather-noomega-nodlogomega-nodeath-all-nobs-nodel-combo.csv',
                  'epimodel-meta-0314-noweather-noomega-nodeath-all-nobs-nodel.csv', ##
                  'epimodel-meta-0314-noweather-all-nobs-nodel.csv',
                  'epimodel-meta-0314-noprior-noomega-all-nobs-nodel-combo.csv', 
                  'epimodel-meta-0314-noprior-all-nobs-nodel.csv', 'epimodel-meta-0314-full3-all-nobs-nodel.csv')
result.names <- c("Kucharski et al.", " + Variable reporting rate", " + Handle deaths",
                  " + Weather effects on transmission", " + Weather effects on detection",
                  " + OLS-based priors")
                  

coeff.names <- list('alpha'='Time Variance (α)', 'invsigma'='Incubation Period (1/σ)', 'invgamma'='Infectious Period (1/γ)', 'omega'='Reporting Rate (ω)', 'error'='Model Error (ε)', 'e.absh'='Beta[Abs. Humid.]', 'e.r'='Beta[Rel. Humid.]', 'e.t2m'='Beta[Surface Temp.]', 'e.tp'='Beta[Total Prec.]', 'o.absh'='Zeta[Abs. Humid.]', 'o.r'='Zeta[Rel. Humid.]', 'o.t2m'='Zeta[Surface Temp.]', 'o.tp'='Zeta[Total Prec.]', 'deathrate'='Death Rate (δ)', 'deathomegaplus'='Death Reporting (λ)', 'deathlearning'='Death rate learning (ν)', 'portion_early'='Portion Early', 'mobility_slope'='Mobility Slope', 'logbeta'='Log Transmission', 'logomega'='Log Reporting', 'eein'='Infected Imports', 'omega0'='Initial Reporting (ω₀)', 'domega'='Reporting Slope', 'invkappa'='Testing Delay (1/κ)', 'invtheta'='Reporting Delay (1/θ)', 'e.ssrd', 'e.ssrd'='Beta[Solar Rad.]', 'e.utci'='Beta[Thermal Index]', 'o.ssrd'='Zeta[Solar Rad.]', 'o.utci'='Zeta[Thermal Index]')
coeff.order <- c('Time Variance (α)', 'Incubation Period (1/σ)', 'Infectious Period (1/γ)', 'Testing Delay (1/κ)', 'Reporting Delay (1/θ)', 'Reporting Rate (ω)', 'Initial Reporting (ω₀)', 'Reporting Slope', 'Death Rate (δ)', 'Death Reporting (λ)', 'Death rate learning (ν)', 'Beta[Abs. Humid.]', 'Beta[Rel. Humid.]', 'Beta[Surface Temp.]', 'Beta[Total Prec.]', 'Beta[Solar Rad.]', 'Beta[Thermal Index]', 'Zeta[Abs. Humid.]', 'Zeta[Rel. Humid.]', 'Zeta[Surface Temp.]', 'Zeta[Total Prec.]', 'Zeta[Solar Rad.]', 'Zeta[Thermal Index]', 'Portion Early', 'Log Transmission', 'Log Reporting', 'Model Error (ε)')

fits <- data.frame()
performs <- data.frame()
for (rr in 1:length(result.files)) {
    filepath <- file.path("../../results", result.files[rr])
    print(filepath)
    df <- read.csv(filepath)
    for (param in unique(df$param)) {
        if (!(coeff.names[[param]] %in% coeff.order))
            next
        
        if ("  " %in% df$regid[df$param == param]) {
            meanmu <- df$mu[df$param == param & df$regid == "  "]
            mean25 <- df$ci25[df$param == param & df$regid == "  "]
            mean75 <- df$ci75[df$param == param & df$regid == "  "]
        } else {
            meanmu <- mean(df$mu[df$param == param], na.rm=T)
            mean25 <- mean(df$ci25[df$param == param], na.rm=T)
            mean75 <- mean(df$ci75[df$param == param], na.rm=T)
        }
        mu25q.raw <- quantile(df$mu[df$param == param & df$group == "Raw"], .25, na.rm=T)
        mu75q.raw <- quantile(df$mu[df$param == param & df$group == "Raw"], .75, na.rm=T)
        mu25q.combined <- quantile(df$mu[df$param == param & df$group == "Combined" & df$regid != "  "], .25, na.rm=T)
        mu75q.combined <- quantile(df$mu[df$param == param & df$group == "Combined" & df$regid != "  "], .75, na.rm=T)
        if ('rhat' %in% names(df))
            rhat <- mean(df$rhat[df$param == param], na.rm=T)
        else
            rhat <- NA

        if (any(result.files == 'epimodel-0803.csv')) {
            if (rr >= which(result.files == 'epimodel-0803.csv') & rr <= which(result.files == 'epimodel-meta-0806-x5.csv')) {
                if (param == 'e.absh')
                    param <- 'e.r'
                else if (param == 'e.r')
                    param <- 'e.t2m'
            }
        }

        ymin <- NA
        ymax <- NA
        if (substring(as.character(param), 1, 2) == 'e.')
            if (!is.na(mean25) && (mean75 > 1 || mean25 < -1 || mu75q.raw > 1 || mu25q.raw < -1)) {
                ymin <- -1
                ymax <- 1
                meanmu <- NA
                mean25 <- NA
                mean75 <- NA
                mu25q <- NA
                mu75q <- NA
            }

        fits <- rbind(fits, data.frame(model=result.names[rr], param, param.name=coeff.names[[param]], meanmu, mean25, mean75, mu25q.raw, mu75q.raw, mu25q.combined, mu75q.combined, ymin, ymax, rhat=rhat))
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
    facet_wrap(~ param.name, scales='free_x', ncol=6) +
    geom_errorbar(aes(ymin=mean25, ymax=mean75)) +
    geom_linerange(aes(ymin=mu25q.combined, ymax=mu75q.combined), col='red', linetype='dashed') +
    # geom_linerange(aes(ymin=ymin, ymax=ymax), col='blue') +
    coord_flip() + scale_y_continuous(expand=c(.01, .01)) + theme_bw() +
    xlab(NULL) + ylab("Parameter Estimate")
ggsave("../../figures/compare-models.pdf", gp, width=10, height=8)

gp <- ggplot(fits[!(fits$param %in% c('logbeta', 'logomega')),], aes(model, meanmu)) +
    facet_wrap(~ param.name, scales='free_x', ncol=6) +
    geom_errorbar(aes(ymin=mu25q.combined, ymax=mu75q.combined)) +
    geom_linerange(aes(ymin=mu25q.raw, ymax=mu75q.raw), col='red', linetype='dashed') +
    # geom_linerange(aes(ymin=ymin, ymax=ymax), col='blue') +
    coord_flip() + scale_y_continuous(expand=c(.01, .01)) + theme_bw() +
    xlab(NULL) + ylab("Parameter Estimate")
ggsave("../../figures/compare-models-0314.pdf", gp, width=13, height=5)

cairo_pdf("../../figures/compare-models-0314.pdf", width=13, height=5)
gp
dev.off()
