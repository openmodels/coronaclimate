library(dplyr)
library(lfe)
library(ggplot2)

df <- read.csv("cases/panel_all.csv")
df$regid <- paste(df$Country, df$Region, df$Locality)
df$days <- as.numeric(difftime(df$Date, "2020-01-01", units='days'))
df$Confirmed[df$Confirmed < .5] <- 0 # stop very low values
df$t2m[df$t2m == 0] <- NA

## Diagnostic
df2 <- df %>% group_by(regid) %>% summarize(index=sum(t2m * tp))
df2$regid[duplicated(df2$index)]

## Construct average over previous period
weather.delay <- c(10, 19)

df$Confirmed.delay <- NA
weathers <- c('absh', 'absh2', 'de', 'de2', 'q', 'q2', 'r', 'r2', 'ssrd', 'ssrd2', 't2m', 't2m2', 'tp', 'tp2')
for (weather in weathers)
    df[, paste0(weather, '.pred')] <- NA

for (regid in unique(df$regid)) {
    print(regid)
    rows <- which(df$regid == regid)
    subdf <- df[rows,]
    if (nrow(subdf) > 200) {
        print(paste("Duplicated region:", regid))
        next
    }
    if (!(subdf$days[1] == 0 && all(round(diff(subdf$days)) == 1))) {
        print(paste("Bad days in region:", regid))
        next
    }

    df$Confirmed.delay[rows] <- c(NA, subdf$Confirmed[-nrow(subdf)])
    for (weather in weathers) {
        values <- convolve(subdf[-nrow(subdf), weather], rep(1, diff(weather.delay)+1) / (diff(weather.delay)+1), type='filter')
        df[rows, paste0(weather, '.pred')] <- c(rep(NA, diff(weather.delay)+1), values[1:(length(values) - diff(weather.delay)-1)], rep(NA, diff(weather.delay)+1))
    }
}

df$dlog <- log(df$Confirmed) - log(df$Confirmed.delay)
df$dlog[!is.finite(df$dlog)] <- NA

## df$tas <- df$t2m.week - 273.15
## df$tas2 <- df$tas^2
## df$tas3 <- df$tas^3

## df$prcp <- df$tp.week
## df$prcp2 <- df$tp.week^2

df$superset <- as.character(df$Country)
df$superset[df$Region == "" & df$Locality == ""] <- "global"

df$week <- floor(df$days / 7)

df$days <- round(df$days, 0)

save(df, file="cases/panel-prepped.RData")
