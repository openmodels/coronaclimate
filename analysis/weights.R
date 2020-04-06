setwd("~/Dropbox/Coronavirus and Climate")

library(PBSmapping)
library(raster)

pop <- raster("~/data/gpw/gpw_v4_global_agg15min_2020_ERA.nc")

pop2 <- as.data.frame(pop, xy=T)
pop2$EID <- 1:nrow(pop2)
names(pop2) <- c("X", "Y", "population", "EID")
events <- as.EventData(pop2, projection=1)

make.weights <- function(shppath, weightpath) {
    shp <- importShapefile(shppath)
    polydata <- attr(shp, 'PolyData')

    mapping <- findPolys(events, shp, maxRows=nrow(pop2))

    ## Construct my own sparse matrix
    indexes <- list()
    weights <- list()

    for (pid in unique(mapping$PID)) {
        indexes[[pid]] <- mapping$EID[mapping$PID == pid]
        weights[[pid]] <- events$population[indexes[[pid]]]
    }

    save(indexes, weights, polydata, file=weightpath)
}

for (adm in 0:2) {
    make.weights(paste0("~/data/gadm/gadm36_levels_simple/adm", adm, ".shp"), paste0("popwt_adm", adm, ".RData"))
}

make.weights("shapefiles/usa/US_county_2000-simple-latlon.shp", "popwt_usa.RData")
