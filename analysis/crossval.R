setwd("~/Dropbox/Coronavirus and Climate")

source("./code/analysis/load.R")

#mod <- felm(dlog ~ tas + tas2 + tas3 + prcp | factor(regid) + factor(regid) : factor(week) + factor(superset) : factor(Date) | 0 | regid, data=df)
#summary(mod)$P.r.squared

## Collect temp. and precip. from last 3 weeks
weathers <- c('t2m', 'tp')
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

## Select 10% of each region with at least 10 values >0
traindata <- c()
testdata <- c()
for (regid in unique(df$regid)) {
    rows <- which(df$regid == regid & !is.na(df$dlog) & !is.na(df$t2m.d21))
    if (sum(df$dlog[rows] > 0, na.rm=T) >= 10) {
        print(regid)
        testdatum <- sample(rows, round(.1 * length(rows)))
        traindata <- c(traindata, rows[!(rows %in% testdatum)])
        testdata <- c(testdata, testdatum)
    }
}

alldata <- c(traindata, testdata)
fl <- list(factor(df$regid[alldata]), factor(paste(df$regid[alldata], df$week[alldata])), factor(paste(df$superset[alldata], df$Date[alldata])))

df$mytas <- df$mytas2 <- df$mytas3 <- df$myprcp <- NA
results <- data.frame()
for (last in 1:21) {
    for (first in last:21) {
        print(c(first, last))
        t2msum <- 0
        tpsum <- 0
        for (delay in last:first) {
            t2msum <- t2msum + df[, paste0('t2m.d', delay)]
            tpsum <- tpsum + df[, paste0('tp.d', delay)]
        }

        df$mytas <- t2msum / (first - last + 1) - 273.15
        df$mytas2 <- df$mytas^2
        df$mytas3 <- df$mytas^3
        df$myprcp <- tpsum / (first - last + 1)

        ## mod <- felm(dlog ~ mytas + mytas2 + mytas3 + myprcp | factor(regid) + factor(regid) : factor(week) + factor(superset) : factor(Date) | 0 | regid, data=df[traindata,])
        ## fes <- getfe(mod)

        dmdf <- demeanlist(df[alldata, c('dlog', 'mytas', 'mytas2', 'mytas3', 'myprcp')], fl)
        mymod <- lm(dlog ~ mytas + mytas2 + mytas3 + myprcp, dmdf[1:length(traindata),])

        preds <- predict(mymod, dmdf[(1+length(traindata)):nrow(dmdf),])

        yytrue <- dmdf$dlog[(1+length(traindata)):nrow(dmdf)]
        rsqr <- sum((preds - mean(yytrue))^2) / sum((yytrue - mean(yytrue))^2)

        results <- rbind(results, data.frame(first, last, rsqr))
    }
}
library(ggplot2)

ggplot(results, aes(first, last, fill=rsqr)) +
    geom_raster() + scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0)) + scale_fill_continuous(name="Crossval\nR-Sqr.", limits=c(0, .004)) +
    theme_bw() + xlab("First month for averaging (days until recorded)") + ylab("Last month for averaging (days until recorded)")

results[which.max(results$rsqr),]
