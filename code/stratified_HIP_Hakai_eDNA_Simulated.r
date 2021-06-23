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

setwd("../")
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

plot(habitat_polygon_features %>% filter(habitat == "bull_kelp"))
plot(habitat_polygon_features %>% filter(habitat == "giant_kelp"))
plot(habitat_polygon_features %>% filter(habitat == "high_rugosity"))
plot(habitat_polygon_features %>% filter(habitat == "giant_kelp"))
plot(habitat_polygon_features %>% filter(habitat == "giant_kelp"))

#read in legacy sites
legacy_sites <- readOGR("./spatial_data/hakai_legacy_sites.gpkg")
legacy_sites <- st_as_sf(legacy_sites)

#stratified BAS sample
#draw 1000 BAS sample points from each layer
N_Zone <- c("high_rugosity" = 20,
            "low_rugosity" = 20,
            "bull_kelp" = 40,
            "giant_kelp" = 60,
            "seagrass" = 40,
            "unclassified" = 20)

sample_frame <- readOGR("./spatial_data/habitat_polygons.gpkg")
sample_frame <- st_as_sf(sample_frame)
pts <- st_cast(sample_frame, "POINT")

# Let's simplify this a little:
set.seed(99)
pts <- pts[sample(nrow(pts), floor(0.1*nrow(pts)), replace = FALSE),]
pts$incl_prob <- 0
pts$trackID <- 1:nrow(pts)

habitats <- names(N_Zone)

pts.list <- list()
for(i in 1:length(habitats)){
  habitat <- habitats[i]
  print(habitat)

  pts.list[[i]] <- filter(pts, habitat == habitats[i])
}


#habitats
# legacy sites to subtract: bull_kelp 5, giant_kelp 7, high_rugosity 0, low_rugosity 0, seagrass 16, unclassified 14 (not in order)
sampleV <- c(13,15,24,46,40,20) # before removing legacy sites it was c(20,20,40,60,40,20) 
selected_pts <- list()

pb <- txtProgressBar(min = 0, max = 1000, style = 3)
for(k in 1:1000)
{
	bb.i <- buildMS(pts, d = 2, FALSE, rotate = TRUE) # build master sample for a sample frame of all pts!!! not pts.i

	for(i in 1:length(habitats)){

	  pts.i <- pts.list[[i]]
	  pts.i <- getIndividualBoxIndices(pts = pts.i, bb = bb.i, size = 100)
	  coords <- st_coordinates(pts.i)
	  ## Bounding box for HIP sample:
	  bb.tmp <- st_bbox(pts.i)
	  bbox <- cbind(c(bb.tmp['xmin'], bb.tmp['ymin']), c(bb.tmp['xmax'], bb.tmp['ymax']))
	  
	  HSS.pts <- getHipSample(X = coords[,1], Y = coords[,2], index = pts.i$HaltonIndex,  
							  N = N_Zone[i], bb = bbox,  base = c(2,3), quiet = TRUE,
							  Ps1 = 0:1, Ps2 = 0:2, hipS1 = sample(0:1,2, replace= FALSE), hipS2 = sample(0:2,3, replace= FALSE))	# Changed the order for fun
	  HSS.pts <- HSS.pts[!is.na(HIPIndex)]
	  n.boxes <- length(table(HSS.pts$HIPIndex))
	  
	  # Chooses a sample from all boxes. Small wrapper function.
	  n <- N_Zone[i]
	  pts.new <- getSamples(HSS.pts, n = n)
	  pts.new$habitat <- habitats[i]
	  tmp.id <- pts.i %>% filter(HaltonIndex %in% pts.new$index) %>% filter(!duplicated(HaltonIndex))
	  pts[pts$trackID %in% tmp.id$trackID, ]$incl_prob <- pts[pts$trackID %in% tmp.id$trackID, ]$incl_prob + 1
	}
   setTxtProgressBar(pb, k)
}

## Do the SRS:
pts$incl_prob_srs <- 0
pb <- txtProgressBar(min = 0, max = 1000, style = 3)
for(k in 1:1000)
{
	for(i in 1:length(habitats))
	{
	  pts.i <- pts.list[[i]]
	  smp <- pts.i[sample(nrow(pts.i), N_Zone[i], replace = FALSE), ]
	  pts[pts$trackID %in% smp$trackID, ]$incl_prob_srs <- pts[pts$trackID %in% smp$trackID, ]$incl_prob_srs + 1  
	}
   setTxtProgressBar(pb, k)
}

par(mfrow = c(1,2))
boxplot(pts[pts$habitat == habitats[1],]$incl_prob/1000, main = "HSS", ylim = c(0,0.03))
abline(h = sum(N_Zone[1])/nrow(pts.list[[1]]), col = "red")
boxplot(pts[pts$habitat == habitats[1],]$incl_prob_srs/k, main = "SRS", ylim = c(0,0.03))
abline(h = sum(N_Zone[1])/nrow(pts.list[[1]]), col = "red")
par(mfrow = c(1,2))
boxplot(pts[pts$habitat == habitats[2],]$incl_prob/1000, main = "HSS", ylim = c(0,0.02))
abline(h = sum(N_Zone[2])/nrow(pts.list[[2]]), col = "red")
boxplot(pts[pts$habitat == habitats[2],]$incl_prob_srs/k, main = "SRS", ylim = c(0,0.02))
abline(h = sum(N_Zone[2])/nrow(pts.list[[2]]), col = "red")
par(mfrow = c(1,2))
boxplot(pts[pts$habitat == habitats[3],]$incl_prob/1000, main = "HSS", ylim = c(0,0.025))
abline(h = sum(N_Zone[3])/nrow(pts.list[[3]]), col = "red")
boxplot(pts[pts$habitat == habitats[3],]$incl_prob_srs/k, main = "SRS", ylim = c(0,0.025))
abline(h = sum(N_Zone[3])/nrow(pts.list[[3]]), col = "red")
par(mfrow = c(1,2))
boxplot(pts[pts$habitat == habitats[4],]$incl_prob/1000, main = "HSS", ylim = c(0,0.04))
abline(h = sum(N_Zone[4])/nrow(pts.list[[4]]), col = "red")
boxplot(pts[pts$habitat == habitats[4],]$incl_prob_srs/k, main = "SRS", ylim = c(0,0.04))
abline(h = sum(N_Zone[4])/nrow(pts.list[[4]]), col = "red")
par(mfrow = c(1,2))
boxplot(pts[pts$habitat == habitats[5],]$incl_prob/1000, main = "HSS", ylim = c(0,0.04))
abline(h = sum(N_Zone[5])/nrow(pts.list[[5]]), col = "red")
boxplot(pts[pts$habitat == habitats[5],]$incl_prob_srs/k, main = "SRS", ylim = c(0,0.04))
abline(h = sum(N_Zone[5])/nrow(pts.list[[5]]), col = "red")
par(mfrow = c(1,2))
boxplot(pts[pts$habitat == habitats[6],]$incl_prob/1000, main = "HSS", ylim = c(0,0.04))
abline(h = sum(N_Zone[6])/nrow(pts.list[[6]]), col = "red")
boxplot(pts[pts$habitat == habitats[6],]$incl_prob_srs/k, main = "SRS", ylim = c(0,0.04))
abline(h = sum(N_Zone[6])/nrow(pts.list[[6]]), col = "red")

