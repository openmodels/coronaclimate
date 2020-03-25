setwd("~/Dropbox/Coronavirus and Climate")

library(PBSmapping)
library(raster)

pop <- raster("~/data/gpw/gpw_v4_global_agg15min_2020_ERA.nc")

pop2 <- as.data.frame(pop, xy=T)
pop2$EID <- 1:nrow(pop2)
names(pop2) <- c("X", "Y", "population", "EID")
events <- as.EventData(pop2, projection=1)

for (adm in 0:2) {
    shp <- importShapefile(paste0("~/data/gadm/gadm36_levels_simple/adm", adm, ".shp"))
    polydata <- attr(shp, 'PolyData')

    mapping <- findPolys(events, shp, maxRows=nrow(pop2))

    ## Construct my own sparse matrix
    indexes <- list()
    weights <- list()

    for (pid in unique(mapping$PID)) {
        indexes[[pid]] <- mapping$EID[mapping$PID == pid]
        weights[[pid]] <- events$population[indexes[[pid]]]
    }

    save(indexes, weights, polydata, file=paste0("popwt_adm", adm, ".RData"))
}
