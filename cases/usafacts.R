setwd("~/Dropbox/Coronavirus and Climate/cases/usafacts")

library(reshape2)
library(dplyr)

confirmed <- read.csv("covid_confirmed_usafacts.csv")
deaths <- read.csv("covid_deaths_usafacts.csv")

## Distribution "unallocated" by population

load("../../weather/popwt_usa.RData")
polydata$fips <- as.numeric(as.character(polydata$NHGISST)) * 100 + as.numeric(as.character(polydata$NHGISCTY)) / 10

redist.unallocs <- function(tbl) {
    for (ii in which(tbl$County.Name == "Statewide Unallocated")) {
        county.rows <- which(tbl$State == tbl$State[ii])
        county.weights <- c()
        for (county.row in county.rows) {
            pids <- polydata$PID[!is.na(polydata$fips) & polydata$fips == tbl$countyFIPS[county.row]]
            if (length(pids) == 0)
                county.weights <- c(county.weights, 0)
            else if (length(pids) == 1)
                county.weights <- c(county.weights, sum(weights[[pids]]))
            else
                print("Multiple PIDs.")
        }
        for (cc in 1:length(county.rows)) {
            tbl[county.rows[cc], 5:ncol(tbl)] <- tbl[county.rows[cc], 5:ncol(tbl)] + tbl[ii, 5:ncol(tbl)] * county.weights[cc] / sum(county.weights)
        }
    }
    subset(tbl, tbl$County.Name != "Statewide Unallocated")
}

confirmed2 <- redist.unallocs(confirmed)
deaths2 <- redist.unallocs(deaths)

confirmed3 <- melt(confirmed2, 1:4, variable.name='xdate', value.name='Confirmed')
deaths3 <- melt(deaths2, 1:4, variable.name='xdate', value.name='Deaths')

df <- confirmed3 %>% left_join(deaths3)
df$date <- as.Date(gsub('\\.', '-', substring(df$xdate, 2, nchar(as.character(df$xdate)))), "%m-%d-%Y")
df$Country <- "USA"
df$Recovered <- NA
df$Source <- "usafacts"

df2 <- df[, c('Country', 'State', 'County.Name', 'date', 'Confirmed', 'Deaths', 'Recovered', 'Source', 'countyFIPS')]
names(df2)[1:4] <- c('Country', 'Region', 'Locality', 'Date')
names(df2)[9] <- "fips"

write.csv(df2, "standardised.csv", row.names=F)

## ## Check the data

## cns <- unique(df$County.Name)
## cns[grep("County", cns)] <- NA
## cns[grep("City", cns)] <- NA
## cns[grep("Paris", cns)] <- NA
## cns[grep("city", cns)] <- NA
## cns[grep("Borough", cns)] <- NA
## unique(cns)

## tbl <- table(confirmed$countyFIPS)
## tbl[tbl > 1]


