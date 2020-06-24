
allrecorded <- read.csv(outfile)

## Show all results
ggplot(subset(allrecorded, param == "omega" & Region == "" & Locality == ""), aes(regid, mu, colour=group)) +
    coord_flip() +
    geom_point() + geom_errorbar(aes(ymin=ci25, ymax=ci75)) +
    theme_bw() + theme(axis.text.y=element_text(size=4)) + ylab(NULL) + xlab(NULL) +
    scale_colour_discrete(name=NULL)

ggplot(subset(allrecorded, country == "Global"), aes(param, mu)) +
    coord_flip() +
    geom_point() + geom_errorbar(aes(ymin=ci25, ymax=ci75)) +
    theme_bw() + ylab("Hyper-paramater value and 75% CI") + xlab(NULL)


## Prepare to map

## Grab ALPHA.3 from df
load("cases/panel-prepped.RData")
df2 <- df[, c('regid', 'ALPHA.3')] %>% group_by(regid) %>% summarize(ALPHA.3=ALPHA.3[1])
allrecorded2 <- subset(allrecorded, Region == "" & Locality == "" & group == "Combined") %>% left_join(df2)

## Map to PIDs
shp <- importShapefile("shapefiles/gadm36_levels_simple/adm0.shp")
polydata <- attr(shp, 'PolyData')

allrecorded3 <- allrecorded2 %>% left_join(polydata[, -1], by=c('ALPHA.3'='GID_0'))

shp2 <- shp %>% left_join(allrecorded3[, c('PID', 'metamu', 'metasd')])

ggplot(shp2, aes(X, Y, fill=metamu, group=paste(PID, SID))) +
    geom_polygon() + scale_y_continuous(name=NULL, limits=c(-60, 85), expand=c(0, 0)) +
    scale_x_continuous(name=NULL, expand=c(0, 0)) + theme_bw() + scale_fill_continuous(name=param)


