setwd("~/Dropbox/Coronavirus and Climate")

load('./cases/panel-prepped.RData')
library(dplyr)
library(tidyr)
library(lfe)
library(ggplot2)

df <- df[df$lowest_level == 1, ]

#mod <- felm(dlog ~ tas + tas2 + tas3 + prcp | factor(regid) + factor(regid) : factor(week) + factor(superset) : factor(Date) | 0 | regid, data=df)
#summary(mod)$P.r.squared

## Collect temp. and precip. from last 3 weeks
weathers <- c('t2m', 'tp', 'absh', 'ssrd', 'de', 'q', 'r')
for (delay in 1:21)
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

    for (delay in 1:21) {
        for (weather in weathers) {
            values <- c(rep(NA, delay), subdf[1:(nrow(subdf) - delay), weather])
            df[rows, paste0(weather, '.d', delay)] <- values
        }
    }
}

#df_orig <- df

ma <- function(x, n = 7){stats::filter(x, rep(1 / n, n), sides = 1)}
df_new <- df %>%
    group_by(regid) %>%
    mutate(dlog_ma = ifelse(!is.na(dlog) & sum(!is.na(dlog))>7, ma(dlog), NA)) %>%
    tidyr::fill(dlog_ma, .direction="updown") %>%
    mutate(dlog_ma = ifelse(!is.na(dlog) & dlog_ma != 0, dlog_ma, NA)) %>%
    mutate(growth = ifelse(dlog_ma <= quantile(dlog_ma, c(0.25), na.rm=T), 'q25',
                    ifelse(dlog_ma <= quantile(dlog_ma, c(0.5), na.rm=T), 'q50',
                    ifelse(dlog_ma <= quantile(dlog_ma, c(0.75), na.rm=T), 'q75',
                    ifelse(dlog_ma <= quantile(dlog_ma, c(1.0), na.rm=T), 'q100', NA))))) %>%
    ungroup()

df_new <- df_new %>%
    group_by(regid) %>%
    mutate(daysince1 = days - days[dlog > 0. & !is.na(dlog)][1]) %>%
    mutate(growth2 = ifelse(daysince1 > 0 & daysince1 <= 14, 'd14', 
                    ifelse(daysince1 > 14 & daysince1 <= 28, 'd28',
                    ifelse(daysince1 > 28, 'dX', NA)))) %>%
    ungroup()

df_new <- df_new %>%
    group_by(regid) %>%
    mutate(month = months(strptime(Date, format='%Y-%m-%d'))) %>%
    mutate(growth3 = ifelse(month == 'March', 'march', 
                    ifelse(month == 'April', 'april', NA))) %>%
    ungroup()

df <- base::data.frame(select(df_new, -c('dlog_ma', 'daysince1', 'month')))

results <- data.frame()
for (growth in c('any', unique(na.omit(df$growth)), unique(na.omit(df$growth2)), unique(na.omit(df$growth3)))) {
    print(growth)
    ## Select 10% of each region with at least 10 values >0
    traindata <- c()
    testdata <- c()
    for (regid in unique(df$regid)) {
        if (growth == 'any') {
            rows <- which(df$regid == regid & !is.na(df$dlog) & !is.na(df$t2m.d21))
        } else if (growth %in% unique(df$growth)) {
            rows <- which(df$regid == regid & !is.na(df$dlog) & !is.na(df$t2m.d21) & df$growth == growth)
        } else if (growth %in% unique(df$growth2)) {
            if (growth == 'd28') {
                rows <- which(df$regid == regid & !is.na(df$dlog) & !is.na(df$t2m.d21) & (df$growth2 == 'd14' | df$growth2 == 'd28'))
            } else {
                rows <- which(df$regid == regid & !is.na(df$dlog) & !is.na(df$t2m.d21) & df$growth2 == growth)
            }
        } else if (growth %in% unique(df$growth3)) {
            rows <- which(df$regid == regid & !is.na(df$dlog) & !is.na(df$t2m.d21) & df$growth3 == growth)
        }
        if (sum(df$dlog[rows] > 0, na.rm=T) >= 10) {
            testdatum <- sample(rows, round(.1 * length(rows)))
            traindata <- c(traindata, rows[!(rows %in% testdatum)])
            testdata <- c(testdata, testdatum)
        }
    }

    alldata <- c(traindata, testdata)
    #print(length(alldata))
    df$mytas <- df$mytas2 <- df$myprcp <- df$myprcp2 <- df$myabsh <- df$myabsh2 <- df$myssrd <- df$myssrd2 <- df$myde <- df$myde2 <- df$myq <- df$myq2 <- df$myr <- df$myr2 <- NA

    for (last in 1:21) {
        for (first in last:21) {
            print(c(first, last))
            t2msum <- 0
            tpsum <- 0
            abshsum <- 0
            ssrdsum <- 0
            desum <- 0
            qsum <- 0
            rsum <- 0
            for (delay in last:first) {
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

            ## mod <- felm(dlog ~ mytas + mytas2 + mytas3 + myprcp | factor(regid) + factor(regid) : factor(week) + factor(superset) : factor(Date) | 0 | regid, data=df[traindata,])
            ## fes <- getfe(mod)

            fl <- list(factor(df$regid[alldata]), factor(paste(df$regid[alldata], df$week[alldata])), factor(paste(df$superset[alldata], df$Date[alldata])))
            dmdf <- demeanlist(df[alldata, c('dlog', 'mytas', 'mytas2', 'myprcp', 'myprcp2', 'myabsh', 'myabsh2', 'myssrd', 'myssrd2')], fl)
            mymod <- lm(dlog ~ mytas + mytas2 + myprcp + myprcp2 + myabsh + myabsh2 + myssrd + myssrd2 , dmdf[1:length(traindata),])

            preds <- predict(mymod, dmdf[(1+length(traindata)):nrow(dmdf),])

            yytrue <- dmdf$dlog[(1+length(traindata)):nrow(dmdf)]
            rsqr <- sum((preds - mean(yytrue))^2) / sum((yytrue - mean(yytrue))^2)
            country <- "allcountries"
            results <- rbind(results, data.frame(first, last, rsqr, country, growth))

            countrylist <- c("United States", "Germany", "United Kingdom")
            for (country in countrylist) {
                print(country)
                countrydata <- which(df$Country == country & df$Locality != "")
                alldata_country <- alldata[(alldata %in% countrydata)]
                traindata_country <- traindata[(traindata %in% countrydata)]
                #print(length(alldata_country))

                fl_country <- list(factor(df$regid[alldata_country]), factor(paste(df$regid[alldata_country], df$week[alldata_country])), factor(paste(df$superset[alldata_country], df$Date[alldata_country])))
                dmdf <- demeanlist(df[alldata_country, c('dlog', 'mytas', 'mytas2', 'myprcp', 'myprcp2', 'myabsh', 'myabsh2', 'myssrd', 'myssrd2', 'myde', 'myde2', 'myq', 'myq2', 'myr', 'myr2')], fl_country)
                mymod <- lm(dlog ~ mytas + mytas2 + myprcp + myprcp2 + myabsh + myabsh2 + myssrd + myssrd2 + myde + myde2 + myq + myq2 + myr + myr2, dmdf[1:length(traindata_country),])

                preds <- predict(mymod, dmdf[(1+length(traindata_country)):nrow(dmdf), ])

                yytrue <- dmdf$dlog[(1+length(traindata_country)):nrow(dmdf)]
                rsqr <- sum((preds - mean(yytrue))^2) / sum((yytrue - mean(yytrue))^2)

                results <- rbind(results, data.frame(first, last, rsqr, country, growth))
            }
        }
    }
}
write.csv(results, 'crossval_results_new.csv')


#library(ggplot2)
#
#ggplot(results, aes(first, last, fill=rsqr)) +
#    geom_raster() + scale_x_continuous(expand=c(0, 0)) +
#    scale_y_continuous(expand=c(0, 0)) + scale_fill_continuous(name="Crossval\nR-Sqr.", limits=c(0, .004)) +
#    theme_bw() + xlab("First month for averaging (days until recorded)") + ylab("Last month for averaging (days until recorded)")
#
#results[which.max(results$rsqr),]
