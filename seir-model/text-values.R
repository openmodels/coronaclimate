setwd("~/research/coronavirus/code/seir-model")

source("../configs.R")

df <- read.csv("../../cases/panel_all.csv")
results <- read.csv("../../results/epimodel-meta-0314-noprior-all-nobs-nodel.csv")

weatherscales <- apply(df[, weather], 2, sd)
units <- c('K', 'mm day^{-1}', 'W m^{-2}', 'point')

issig.95 <- function(ci025, ci975)
    sign(ci975) == sign(ci025)

report.95 <- function(mu, ci025, ci975, unit) {
    digits <- max(2, ceiling(log10(abs(mu)) - log10(ci975 - ci025)) + 1)
    paste0("\\num{", format(mu, digits=digits), "} [\\numrange{", format(ci025, digits=digits),
           "}{", format(ci975, digits=digits), "}] \\unit{", unit, "}")
}

for (ww in weather) {
    rows <- subset(results, Country == "" & Region == "" & param %in% c(paste0('e.', ww), paste0('o.', ww)))

    for (rr in 1:nrow(rows)) {
        ## beta f / sd * (sd / N K) * (100 % / f)
        scale <- 100 / weatherscales[weather == ww]
        if (ww == 'tp')
            scale <- scale / (24 * 1000)
        str <- report.95(rows$mu[rr] * scale, rows$ci2.5[rr] * scale,
                         rows$ci97.5[rr] * scale, "\\%")
        print(paste0(rows$param[rr], ": ", ifelse(issig.95(rows$ci2.5[rr], rows$ci97.5[rr]), '***', '')))
        cat(paste0(str, "\n"))
    }
}

report.range <- function(low, high, unit, digits=NA) {
    if (is.na(digits))
        digits <- max(2, ceiling(log10(abs((low + high)/2)) - log10(high - low)) + 1)
    paste0("\\qtyrange{", format(low, digits=digits), "}{", format(high, digits=digits), "}{", unit, "}")
}

for (ww in weather) {
    for (prefix in c('e.', 'o.')) {
        rows <- subset(results, Country != "" & Region == "" & Locality == "" & group == "Combined" &
                                param == paste0(prefix, ww))

        ## beta f / sd * (sd / N K) * (100 % / f)
        scale <- 100 / weatherscales[weather == ww]
        if (ww == 'tp')
            scale <- scale / (24 * 1000)
        str <- report.range(quantile(rows$mu, .025), quantile(rows$mu, .975), "\\%")
        print(paste0(prefix, ww))
        cat(paste0(str, "\n"))
    }
}

## Model 3 months at +1 C
## FIRST: Construct global model

source("forward-0314.R")

rows <- subset(results, Country == "" & Region == "" & Locality == "" & group == "Combined")
