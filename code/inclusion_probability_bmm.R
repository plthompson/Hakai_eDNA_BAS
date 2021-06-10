rm(list = ls())

#install.packages('bcmapsdata', repos='https://bcgov.github.io/drat/')

library(BASMasterSample)
#library(bcmaps)
#library(bcmapsdata)
library(sf)
library(rgdal)
library(here)
library(tidyverse)

####eDNA sites####

sites <- readOGR("./Hakai_eDNA_BAS/sites0603.shp")
sites <- st_as_sf(sites, wtk = geometry)

unclassified.sites <- filter(sites, habitat == "unclassified")

research_area <- readOGR("./Hakai_eDNA_BAS/spatial_data/calvert_research_area.gpkg")
research_area <- st_as_sf(research_area, wtk = geometry)

habitat_line_features <- readOGR("./Hakai_eDNA_BAS/spatial_data/habitat_line_features.gpkg")
habitat_line_features <- st_as_sf(habitat_line_features)

habitat_polygons <- readOGR("./Hakai_eDNA_BAS/spatial_data/habitat_polygons.gpkg")
habitat_polygons <- st_as_sf(habitat_polygons)

n.smp <- 254
NRun <- 100
incl.count <- data.frame(name = row.names(sites), count = 0)
for(i in 1:NRun)
{
  bb.sites <- buildMS(habitat_polygons, d = 2, FALSE, rotate = TRUE) # build master sample
  pts.hi <- getIndividualBoxIndices(pts = sites, bb = bb.sites, size = 100) # get HIs for each point
  t1 <- masterSample(habitat_polygons, bb = bb.sites) # draw new set of points
  pts <- getIndividualBoxIndices(pts = t1, bb = bb.sites, size = 100) #get HIs for new set of points
  incl.count$count <- incl.count$count + ifelse(pts.hi$HaltonIndex %in% pts$HaltonIndex, 1,0) # tally counts for each site by HI
} 
 
  
#sum(incl.count$count)


boxplot(incl.count[,2]/NRun, main = "HaltonIndex w/ random Rotation", ylim = c(0.1, 1))
boxplot(incl.count[,3]/NRun, main = "Simple Random Sample", ylim = c(0.1, 1)) 



####cities example####

cities <- get_layer("bc_cities")

n.smp <- 20
NRun <- 1000
incl.count <- data.frame(name = cities$NAME, count = 0, srs = 0)
for(i in 1:NRun)
{
  bb.city <- buildMS(cities, d = 2, FALSE)
  pts.hi <- getIndividualBoxIndices(pts = cities, bb = bb.city)
  smp.indx <- which(rank(pts.hi$HaltonIndex) <= n.smp)
  incl.count$count[smp.indx] <- incl.count$count[smp.indx] + 1
  srs.i <- sample(nrow(incl.count), n.smp, replace = FALSE)
  incl.count$srs[srs.i] <- incl.count$srs[srs.i] + 1
}

par(mfrow = c(1,2))
boxplot(incl.count[,2]/NRun, main = "HaltonIndex w/ random Rotation", ylim = c(0.05, 0.2))
boxplot(incl.count[,3]/NRun, main = "Simple Random Sample", ylim = c(0.05, 0.2))


