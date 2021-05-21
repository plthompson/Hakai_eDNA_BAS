# code by Patrick Thompson
# April 2021
# patrick.thompson@dfo-mpo.gc.ca

#load packages#####
library(BASMasterSample)
library(rgdal)
library(sf)
library(tidyverse)
library(data.table)
library(raster)

#change plot theme####
source("./code/functions/plot_theme.R")
theme_set(plt_theme)

#load research area####
research_area <- readOGR("./spatial_data/calvert_research_area.gpkg")
research_area <- st_as_sf(research_area, wtk = geometry)
area_plot <- ggplot()+
  geom_sf(data = research_area, fill = NA)
#area_plot

#read in habitat lines
habitat_line_features <- readOGR("./spatial_data/habitat_line_features.gpkg")
habitat_line_features <- st_as_sf(habitat_line_features)


#stratified BAS sample
#draw 1000 BAS sample points from each layer
N_Zone <- c("high_rugosity" = 40,
            "low_rugosity" = 20,
            "bull_kelp" = 20,
            "giant_kelp" = 20,
            "seagrass" = 40,
            "unclassified" = 60)


#create Halton boxes#####
# note - bounding box should be all of BC to make boxes align with greater BC BAS
# however this is slow and so for now we will use the study region as the bounding box
#make bounding box
# bb <- buildMS(research_area, d = 2, showOutput = FALSE, rotate = FALSE)
# MARINE MS, need to add the seed to it for it to work.
bb <- getBB()
attr(bb, "seed") <- getSeed()

box_size <- 100 #chose size of halton box

halton_boxes <- point2Frame(pts = habitat_line_features, bb = bb, size = box_size)

#area_plot+
#  geom_sf(data = habitat_line_features, size = 0.3) +
#  geom_sf(data = halton_boxes, fill = NA, size = 0.3, color = 1)

#perform Halton Iterative Partitioning####
#library(devtools)
#install_github("paul-vdb/NWTMasterSample", build_vignettes = TRUE)
source('./code/functions/HSS.r')

coords <- st_coordinates(st_centroid(halton_boxes))
bb.tmp <- st_bbox(bb)
bbox <- cbind(c(bb.tmp['xmin'], bb.tmp['ymin']), c(bb.tmp['xmax'], bb.tmp['ymax']))

# Paul: Let's add more boxes to show better spatial balance.
HSS.pts <- getHipSample(X = coords[,1], Y = coords[,2], index = halton_boxes$HaltonIndex,
                        N = 70, bb = bbox,  base = c(2,3), quiet = TRUE,
                        Ps1 = 0:1, Ps2 = 0:2, hipS1 = 0:1, hipS2 = c(0,2,1))	# Changed the order for fun
n.boxes <- length(table(HSS.pts$HIPIndex))

# Chooses a sample from all boxes. Small wrapper function.
n <- 20
pts.new <- getSamples(HSS.pts, n = n)
# Bounding box to clip the HIP boxes to.
bb.tmp <- st_bbox(research_area)
bbox2 <- cbind(c(bb.tmp['xmin'], bb.tmp['ymin']), c(bb.tmp['xmax'], bb.tmp['ymax']))
HIPBoxes <- st_as_sf(getHIPBoxes(hip = HSS.pts, bb = bbox2, n = n.boxes, base = c(2,3)))
HIPBoxes <- st_set_crs(HIPBoxes, value = st_crs(halton_boxes)) #something is wrong here since the HIP boxes don't line up properly

halton_boxes_select <- halton_boxes %>% filter(HaltonIndex %in% pts.new$index)

area_plot+
  geom_sf(data = habitat_line_features, size = 0.3) +
  geom_sf(data = HIPBoxes, fill = NA)+
  geom_sf(data = halton_boxes_select, fill = NA, size = 1, color = 2)

# Check Spatial Balance:
library(rgeos)
sp.var <- sp.var.srs <- NULL

# HSS Sample with varying numbers of partitions but no other randomness.
for(B in c(1, 6, 12, 24, 32, 72, 144)){
	HSS.pts <- getHipSample(X = coords[,1], Y = coords[,2], index = halton_boxes$HaltonIndex,
                        N = B, bb = bbox,  base = c(2,3), quiet = TRUE,
                        Ps1 = 0:1, Ps2 = 0:2, hipS1 = 0:1, hipS2 = c(0,2,1))	# Changed the order for fun
	sp.var <- c(sp.var, getBalance(HSS.pts, n = boxes))
}

# Simple Random Sample
for(i in 1:100){
	smp <- HSS.pts[sample(nrow(HSS.pts), boxes)]
	smp <- SpatialPoints(cbind(smp$X, smp$Y))
	sp.var.srs <- c(sp.var.srs, getBalance(HSS.pts, n = boxes, smp = smp))
}

# Plot of the spatial balance.
plot(c(1, 6, 12, 24, 32, 72, 144), sp.var, type = 'l', ylim = c(0, 60), main = "Variance of Voronoi Polygons")
abline(h = mean(sp.var.srs), col = "red")	# Simple Random Sample
