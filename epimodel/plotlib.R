labelmap <- list('mobility_slope'="Mobility Adjustment",
                 'alpha'="Gradual Adjustment Rate",
                 'invsigma'="Incubation Period (days)",
                 'invgamma'="Infectious Period (days)",
                 'invkappa'="Pre-Testing Period (days)",
                 'invtheta'="Reporting Delay (days)",
                 'logbeta'="Log Transmission Rate",
                 'omega'="Realised Reporting Rate",
                 'deathrate'="Death Rate",
                 'deathomegaplus'="Extra Record of Deaths",
                 'deathlearning'="Death Learning Rate",
                 'portion_early'="Portion Reported Early",
                 'e.t2m'="Air Temperature Trans.",
                 'e.tp'="Total Precipitation Trans.",
                 'e.r'="Relative Humidity Trans.",
                 'e.absh'="Absolute Humidity Trans.",
                 'e.ssrd'="Solar Radiation Trans.",
                 'e.utci'="Thermal Discomfort Trans.",
                 'o.t2m'="Air Temperature Detect",
                 'o.tp'="Total Precipitation Detect",
                 'o.r'="Relative Humidity Detect",
                 'o.absh'="Absolute Humidity Detect",
                 'o.ssrd'="Solar Radiation Detect",
                 'o.utci'="Thermal Discomfort Detect",
                 'error'="Model Error",
                 'logomega'="Log Reporting Rate",
                 'eein'="Exposed Imports",
                 'rhat'="Bayesian convergence")
paramorder <- c('alpha', 'invgamma', 'invsigma', 'invkappa', 'invtheta',
                'mobility_slope', 'omega', 'portion_early',
                'deathrate', 'deathomegaplus', 'deathlearning',
                'logbeta', 'logomega', 'eein',
                'e.absh', 'e.r', 'e.t2m', 'e.tp', 'e.ssrd', 'e.utci',
                'o.absh', 'o.r', 'o.t2m', 'o.tp', 'o.ssrd', 'o.utci', 'error', 'rhat')

## Prefer to use this, since we've messed up with factor params before
get.param.labels <- function(params) {
    param.labels <- sapply(params, function(param) labelmap[[as.character(param)]])
    factor(param.labels, levels=rev(sapply(paramorder, function(param) labelmap[[param]])))
}
