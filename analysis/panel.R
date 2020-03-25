setwd("~/Dropbox/Coronavirus and Climate")

library(ncdf4)
library(PBSmapping)

df <- read.csv("cases/john-hopkins/panel_john-hopkins.csv")
df$regid <- paste(df$Country, df$Region, df$Locality)

nc <- nc_open("~/data/ecmwf-covid/era5_2020-01-01_2020-03-17_t2m-tp-d2m-sp-ssrd_daymean.nc")
## ignore top and bottom lat of nc
var.t2m <- ncvar_get(nc, 't2m')[,2:720,]
var.tp <- ncvar_get(nc, 'tp')[,2:720,]
var.time <- ncvar_get(nc, 'time')

load("popwt_adm0.RData")

weather <- data.frame()
for (regid in unique(df$regid)) {
    print(regid)
    rows <- which(df$regid == regid)
    rr <- rows[1]
    if (df$Region[rr] == '' || is.na(df$Region)) {
        country <- df$Country[rr]
        if (country == 'US')
            country <- "United States"
        adm <- 0
        pid <- polydata$PID[polydata$NAME_0 == as.character(df$Country[rr])]
        if (length(pid) == 0) {
            print(paste("Cannot find country", df$Country[rr]))
            next
        }
    } else {
        print(paste("Cannot process region", regid))
        next
    }

    for (tt in 1:length(var.time)) {
        slice <- t(var.t2m[,, tt])
        avg.t2m <- sum(slice[indexes[[pid]]] * weights[[pid]], na.rm=T) / sum(weights[[pid]], na.rm=T)

        slice <- t(var.tp[,, tt])
        avg.tp <- sum(slice[indexes[[pid]]] * weights[[pid]], na.rm=T) / sum(weights[[pid]], na.rm=T)

        date <- substring(as.POSIXct(var.time[tt]*3600,origin='1900-01-01 00:00'), 1, 10)
        weather <- rbind(weather, data.frame(date, adm, pid, regid, t2m=avg.t2m, tp=avg.tp))
    }
}

write.csv(weather, "~/data/ecmwf-covid/weather.csv", row.names=F)


