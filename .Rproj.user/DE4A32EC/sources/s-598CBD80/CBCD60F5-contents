# code by Patrick Thompson (modified by Ben Millard-Martin)
# April 2021
# patrick.thompson@dfo-mpo.gc.ca

rm(list = ls())

#load packages#####
library(BASMasterSample)
library(rgdal)
library(sf)
library(tidyverse)
library(data.table)
library(raster)
library(cowplot)
library(here)

source("./code/functions/HSS.r")

#change plot theme####
source("./code/functions/plot_theme.R")
theme_set(plt_theme)

#load research area####
research_area <- readOGR("./spatial_data/calvert_research_area.gpkg")
research_area <- st_as_sf(research_area)

#read in habitat polygons
habitat_polygon_features <- readOGR("./spatial_data/habitat_discrete_legacyBuff.gpkg")
habitat_polygon_features <- st_as_sf(habitat_polygon_features)

#read in legacy sites
legacy_sites <- readOGR("./spatial_data/hakai_legacy_sites.gpkg")
legacy_sites <- st_as_sf(legacy_sites)

#stratified BAS sample
#draw 1000 BAS sample points from each layer
N_Zone <- c("high_rugosity" = 1000000,
            "low_rugosity" = 1000000,
            "bull_kelp" = 1000000,
            "giant_kelp" = 1000000,
            "seagrass" = 1000000,
            "unclassified" = 1000000)

areaBAS <- masterSample(shp = habitat_polygon_features, N = N_Zone, stratum = "habitat")
areaBAS$layer <- c(rep(names(N_Zone)[1], N_Zone[1]),
                   rep(names(N_Zone)[2], N_Zone[2]),
                   rep(names(N_Zone)[3], N_Zone[3]),
                   rep(names(N_Zone)[4], N_Zone[4]),
                   rep(names(N_Zone)[5], N_Zone[5]),
                   rep(names(N_Zone)[6], N_Zone[6]))

#double check that this order doesn't change
areaBAS <- mutate(areaBAS, habitat = recode(areaBAS$layer,
                           "giant_kelp" = 1, "bull_kelp"=2,"seagrass"=3,"unclassified"=4,"high_rugosity"=5,"low_rugosity"=6))

#get bounding box for BC                                        
# MARINE MS, need to add the seed to it for it to work.           ###not quite sure what this means
bb <- getBB()
attr(bb, "seed") <- getSeed()


habitats <- unique(habitat_polygon_features$habitat)
#habitats
# legacy sites to subtract: bull_kelp 5, giant_kelp 7, high_rugosity 0, low_rugosity 0, seagrass 16, unclassified 14 (not in order)
sampleV <- c(13,15,24,46,40,20) # before removing legacy sites it was c(20,20,40,60,40,20) 
selected_pts <- list()

for(i in 1:length(habitats)){
  habitat <- habitats[i]
  print(habitat)

  hab_pts <- filter(areaBAS, habitat == i)
  coords <- st_coordinates(hab_pts)
  bb.tmp <- st_bbox(bb)
  bbox <- cbind(c(bb.tmp['xmin'], bb.tmp['ymin']), c(bb.tmp['xmax'], bb.tmp['ymax']))
  
  HSS.pts <- getHipSample(X = coords[,1], Y = coords[,2], index = hab_pts$SiteID,
                          N = sampleV[i]*2, bb = bbox,  base = c(2,3), quiet = TRUE,
                          Ps1 = 0:1, Ps2 = 0:2, hipS1 = 0:1, hipS2 = c(0,2,1))	# Changed the order for fun
  n.boxes <- length(table(HSS.pts$HIPIndex))
  
  # Chooses a sample from all boxes. Small wrapper function.
  n <- sampleV[i]
  pts.new <- getSamples(HSS.pts, n = n)
  pts.new$habitat <- habitats[i]
  
  selected_pts[[i]] <- pts.new 
}

selected_pts <- rbind(selected_pts[[1]],selected_pts[[2]], selected_pts[[3]], selected_pts[[4]], selected_pts[[5]], selected_pts[[6]])

selected_points = st_as_sf(selected_pts, coords = c("X", "Y"), 
                 crs = 3005, agr = "constant")


selected_points$type <- "BAS"
legacy_sites$type <- "legacy"

t1 <- selected_points[,c(13,14,15)]
t2 <- legacy_sites[,c(2,3,4)]
t3 <- rbind(t1, t2)

ggplot() +
  geom_sf(data = research_area, fill = NA)+
  geom_sf(data = habitat_polygon_features, size = 0.3) +
  geom_sf(data = selected_points, aes(fill = habitat), color = 1, size = 3, pch = 21)+
  facet_wrap(~habitat)
ggsave("./figures/HIPpoints1000000_facet.pdf", height = 8, width = 8)

ggplot() +
  geom_sf(data = research_area, fill = NA)+
  geom_sf(data = habitat_polygon_features, size = 0.3) +
  geom_sf(data = t3, aes(color = habitat, shape = type), size = 2)+
  theme(legend.position = c(1,1), legend.justification = c(1,1)) 
ggsave("./figures/HIPpoints1000000.pdf", height = 8, width = 8)




