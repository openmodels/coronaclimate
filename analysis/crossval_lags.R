setwd("~/Dropbox/Coronavirus and Climate")

#source("./code/analysis/load.R")
load('./cases/panel-prepped.RData')
library(dplyr)
library(tidyr)
library(lfe)
library(ggplot2)

#mod <- felm(dlog ~ tas + tas2 + tas3 + prcp | factor(regid) + factor(regid) : factor(week) + factor(superset) : factor(Date) | 0 | regid, data=df)
#summary(mod)$P.r.squared

df <- df[df$lowest_level == 1, ]

## Collect temp. and precip. from last 3 weeks
weathers <- c('t2m', 'tp', 'absh', 'ssrd', 'de', 'q', 'r')
for (delay in 1:30)
    for (weather in weathers)
        df[, paste0(weather, '.d', delay)] <- NA

for (regid in unique(df$regid)) {
    print(regid)
    rows <- which(df$regid == regid)
    subdf <- df[rows,]
    if (nrow(subdf) > 130) {
        print(paste("Duplicated region:", regid))
        next
    }
    stopifnot(subdf$days[1] == 0 && all(diff(subdf$days) == 1))

    for (delay in 1:30) {
        for (weather in weathers) {
            values <- c(rep(NA, delay), subdf[1:(nrow(subdf) - delay), weather])
            df[rows, paste0(weather, '.d', delay)] <- values
        }
    }
}
save(df, file='~/temp.RData')

traindata <- c()
testdata <- c()
for (regid in unique(df$regid)) {
    rows <- which(df$regid == regid & !is.na(df$dlog) & !is.na(df$t2m.d21))
    if (sum(df$dlog[rows] > 0, na.rm=T) >= 10) {
        testdatum <- sample(rows, round(.1 * length(rows)))
        traindata <- c(traindata, rows[!(rows %in% testdatum)])
        testdata <- c(testdata, testdatum)
    }
}

lagfile <- './cases/onset-symptoms_to_confirmation.csv'
df.lags <- read.csv(lagfile)
nobs.min <- 2
lag.max <- 10

df.lags <- df.lags[df.lags$median <= lag.max, ]
countrylist <- as.character(unique(df$Country)[unique(df$Country) %in% df.lags$X[df.lags$count >= nobs.min]])
rows.restofworld <- which(!df$Country %in% countrylist)

results <- data.frame()
alldata <- c(traindata, testdata)
#print(length(alldata))
df$mytas <- df$mytas2 <- df$myprcp <- df$myprcp2 <- df$myabsh <- df$myabsh2 <- df$myssrd <- df$myssrd2 <- df$myde <- df$myde2 <- df$myq <- df$myq2 <- df$myr <- df$myr2 <- NA

for (lag.mode in c('different', 'same')) {
    print(lag.mode)
    for (last in 1:20) {
        for (first in last:20) {
            print(c(first, last))
            if (lag.mode == 'same') {
                lag.reporting <- df.lags[df.lags$X == 'global_mean-of-countries', 'median']
                t2msum <- 0
                tpsum <- 0
                abshsum <- 0
                ssrdsum <- 0
                desum <- 0
                qsum <- 0
                rsum <- 0
                for (delay in last:(first+lag.reporting)) {
                    t2msum <- t2msum + df[, paste0('t2m.d', delay)]
                    tpsum <- tpsum + df[, paste0('tp.d', delay)]
                    abshsum <- abshsum + df[, paste0('absh.d', delay)]
                    ssrdsum <- ssrdsum + df[, paste0('ssrd.d', delay)]
                    desum <- desum + df[, paste0('de.d', delay)]
                    qsum <- qsum + df[, paste0('q.d', delay)]
                    rsum <- rsum + df[, paste0('r.d', delay)]
                }
                df$mytas <- t2msum / (first - last + 1) - 273.15
                df$mytas2 <- df$mytas^2
                df$myprcp <- tpsum / (first - last + 1)
                df$myprcp2 <- df$myprcp^2
                df$myabsh <- abshsum / (first - last + 1)
                df$myabsh2 <- df$myabsh^2
                df$myssrd <- ssrdsum / (first - last + 1)
                df$myssrd2 <- df$myssrd^2
                df$myde <- desum / (first - last + 1)
                df$myde2 <- df$myde^2
                df$myq <- qsum / (first - last + 1)
                df$myq2 <- df$myq^2
                df$myr <- rsum / (first - last + 1)
                df$myr2 <- df$myr^2
            } else {
                for (country in c('restofworld', countrylist)) {
                    print(country)
                    if (country == 'restofworld') {
                        lag.reporting <- df.lags[df.lags$X == 'global_mean-of-countries', 'median']
                        rows <- rows.restofworld
                    } else {
                        lag.reporting <- df.lags[df.lags$X == country, 'median']
                        rows <- which(df$Country == country)
                    }
                    print(lag.reporting)
                    t2msum <- 0
                    tpsum <- 0
                    abshsum <- 0
                    ssrdsum <- 0
                    desum <- 0
                    qsum <- 0
                    rsum <- 0
                    for (delay in last:(first+lag.reporting)) {
                        t2msum <- t2msum + df[rows, paste0('t2m.d', delay)]
                        tpsum <- tpsum + df[rows, paste0('tp.d', delay)]
                        abshsum <- abshsum + df[rows, paste0('absh.d', delay)]
                        ssrdsum <- ssrdsum + df[rows, paste0('ssrd.d', delay)]
                        desum <- desum + df[rows, paste0('de.d', delay)]
                        qsum <- qsum + df[rows, paste0('q.d', delay)]
                        rsum <- rsum + df[rows, paste0('r.d', delay)]
                    }
                    df[rows, ]$mytas <- t2msum / (first - last + 1) - 273.15
                    df[rows, ]$mytas2 <- df[rows, ]$mytas^2
                    df[rows, ]$myprcp <- tpsum / (first - last + 1)
                    df[rows, ]$myprcp2 <- df[rows, ]$myprcp^2
                    df[rows, ]$myabsh <- abshsum / (first - last + 1)
                    df[rows, ]$myabsh2 <- df[rows, ]$myabsh^2
                    df[rows, ]$myssrd <- ssrdsum / (first - last + 1)
                    df[rows, ]$myssrd2 <- df[rows, ]$myssrd^2
                    df[rows, ]$myde <- desum / (first - last + 1)
                    df[rows, ]$myde2 <- df[rows, ]$myde^2
                    df[rows, ]$myq <- qsum / (first - last + 1)
                    df[rows, ]$myq2 <- df[rows, ]$myq^2
                    df[rows, ]$myr <- rsum / (first - last + 1)
                    df[rows, ]$myr2 <- df[rows, ]$myr^2
                }
            }

            ## mod <- felm(dlog ~ mytas + mytas2 + mytas3 + myprcp | factor(regid) + factor(regid) : factor(week) + factor(superset) : factor(Date) | 0 | regid, data=df[traindata,])
            ## fes <- getfe(mod)

            fl <- list(factor(df$regid[alldata]), factor(paste(df$regid[alldata], df$week[alldata])), factor(paste(df$superset[alldata], df$Date[alldata])))
            dmdf <- demeanlist(df[alldata, c('dlog', 'mytas', 'mytas2', 'myprcp', 'myprcp2', 'myabsh', 'myabsh2', 'myssrd', 'myssrd2', 'myde', 'myde2', 'myq', 'myq2', 'myr', 'myr2')], fl)
            mymod <- lm(dlog ~ mytas + mytas2 + myprcp + myprcp2 + myabsh + myabsh2 + myssrd + myssrd2 + myde + myde2 + myq + myq2 + myr + myr2, dmdf[1:length(traindata),])

            preds <- predict(mymod, dmdf[(1+length(traindata)):nrow(dmdf),])

            yytrue <- dmdf$dlog[(1+length(traindata)):nrow(dmdf)]
            rsqr <- sum((preds - mean(yytrue))^2) / sum((yytrue - mean(yytrue))^2)
            results <- rbind(results, data.frame(first, last, rsqr, lag.mode, nobs.min))
        }
    }
}
write.csv(results, 'crossval_results_lags3.csv')

#library(ggplot2)
#
#ggplot(results, aes(first, last, fill=rsqr)) +
#    geom_raster() + scale_x_continuous(expand=c(0, 0)) +
#    scale_y_continuous(expand=c(0, 0)) + scale_fill_continuous(name="Crossval\nR-Sqr.", limits=c(0, .004)) +
#    theme_bw() + xlab("First month for averaging (days until recorded)") + ylab("Last month for averaging (days until recorded)")
#
#results[which.max(results$rsqr),]
