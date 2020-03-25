setwd("~/Dropbox/Coronavirus and Climate")

library(dplyr)
library(lfe)
library(ggplot2)

df <- read.csv("cases/john-hopkins/panel_john-hopkins.csv")
df$regid <- paste(df$Country, df$Region, df$Locality)

weather <- read.csv("~/data/ecmwf-covid/weather.csv")

## Construct average over previous week
smw <- data.frame()
dates <- unique(weather$date)
for (regid in unique(weather$regid)) {
    print(regid)
    subw <- weather[weather$regid == regid,]
    for (tt in 8:length(dates)) {
        t2m <- mean(subw$t2m[(tt-7):(tt - 1)])
        tp <- sum(subw$tp[(tt-7):(tt - 1)])

        smw <- rbind(smw, data.frame(regid, date=dates[tt], t2m, tp))
    }
}
smw$date <- as.character(smw$date)

df$Date.delay <- as.character(as.Date(df$Date) - 1)

df2 <- df %>% left_join(df, by=c('regid', 'Date.delay'='Date'), suffix=c('', '.delay')) %>%
    left_join(smw, by=c('regid', 'Date'='date'))

df2$dlog <- log(df2$Confirmed) - log(df2$Confirmed.delay)
df2$dlog[!is.finite(df2$dlog)] <- NA

df2$tas <- df2$t2m - 273.15
df2$tas2 <- df2$tas^2
df2$tp2 <- df2$tp^2

mod <- felm(dlog ~ tas + tas2 + tp + tp2 | regid + Date | 0 | regid, data=df2)
mod <- felm(dlog ~ tas + tas2 + tp | regid + Date | 0 | regid, data=df2)

source("code/analysis/lib.R")

preddf <- data.frame(tmean=seq(-10, 40), tp=0, tp2=0)
preddf$tas <- preddf$tmean - 20
preddf$tas2 <- preddf$tmean^2 - 20^2

plot.doseresp(preddf, "Average temperature", "Change in growth rate", indep='tmean')
