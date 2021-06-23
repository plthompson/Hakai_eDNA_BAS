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

plot(habitat_polygon_features)

#read in legacy sites
legacy_sites <- readOGR("./spatial_data/hakai_legacy_sites.gpkg")
legacy_sites <- st_as_sf(legacy_sites)

#stratified BAS sample
#draw 1000 BAS sample points from each layer
N_Zone <- c("high rugosity" = 20,
            "low rugosity" = 20,
            "biogenic" = 140,
            "unclassified" = 20)

load("./spatial_data/box_centroids.RData")
# Sample frame is all_centroids

# Let's simplify this a little:
# set.seed(99)
# pts <- pts[sample(nrow(pts), floor(0.1*nrow(pts)), replace = FALSE),]
all_centroids$incl_prob <- 0
all_centroids$trackID <- 1:nrow(all_centroids)

habitats <- names(N_Zone)

pts.list <- list()
for(i in 1:length(habitats)){
  habitat <- habitats[i]
  print(habitat)

  pts.list[[i]] <- filter(all_centroids, habitat == habitats[i])
}

iters <- 1000
pb <- txtProgressBar(min = 0, max = iters, style = 3)
for(k in 1:iters)
{
	bb.i <- buildMS(all_centroids, d = 2, FALSE, rotate = TRUE) # build master sample for a sample frame of all pts!!! not pts.i

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
	  all_centroids[all_centroids$trackID %in% tmp.id$trackID, ]$incl_prob <- all_centroids[all_centroids$trackID %in% tmp.id$trackID, ]$incl_prob + 1
	}
   setTxtProgressBar(pb, k)
}

## Do the SRS:
all_centroids$incl_prob_srs <- 0
pb <- txtProgressBar(min = 0, max = 1000, style = 3)
for(k in 1:iters)
{
	for(i in 1:length(habitats))
	{
	  pts.i <- pts.list[[i]]
	  smp <- pts.i[ceiling(runif(N_Zone[i], 0, nrow(pts.i))), ]
	  all_centroids[all_centroids$trackID %in% smp$trackID, ]$incl_prob_srs <- all_centroids[all_centroids$trackID %in% smp$trackID, ]$incl_prob_srs + 1  
	}
   setTxtProgressBar(pb, k)
}

par(mfrow = c(4,2))
boxplot(all_centroids[all_centroids$habitat == habitats[1],]$incl_prob/iters, main = "HSS", ylim = c(0,0.01), ylab = "Inclusion Probability", xlab = paste(habitats[1]))
abline(h = sum(N_Zone[1])/nrow(pts.list[[1]]), col = "red")
boxplot(all_centroids[all_centroids$habitat == habitats[1],]$incl_prob_srs/iters, main = "SRS", ylim = c(0,0.01), ylab = "Inclusion Probability", xlab = paste(habitats[1]))
abline(h = sum(N_Zone[1])/nrow(pts.list[[1]]), col = "red")
# par(mfrow = c(1,2))
boxplot(all_centroids[all_centroids$habitat == habitats[2],]$incl_prob/iters, main = "HSS", ylim = c(0,0.01), ylab = "Inclusion Probability", xlab = paste(habitats[2]))
abline(h = sum(N_Zone[2])/nrow(pts.list[[2]]), col = "red")
boxplot(all_centroids[all_centroids$habitat == habitats[2],]$incl_prob_srs/iters, main = "SRS", ylim = c(0,0.01), ylab = "Inclusion Probability", xlab = paste(habitats[2]))
abline(h = sum(N_Zone[2])/nrow(pts.list[[2]]), col = "red")
# par(mfrow = c(1,2))
boxplot(all_centroids[all_centroids$habitat == habitats[3],]$incl_prob/iters, main = "HSS", ylim = c(0,0.025), ylab = "Inclusion Probability", xlab = paste(habitats[3]))
abline(h = sum(N_Zone[3])/nrow(pts.list[[3]]), col = "red")
boxplot(all_centroids[all_centroids$habitat == habitats[3],]$incl_prob_srs/iters, main = "SRS", ylim = c(0,0.025), ylab = "Inclusion Probability", xlab = paste(habitats[3]))
abline(h = sum(N_Zone[3])/nrow(pts.list[[3]]), col = "red")
# par(mfrow = c(1,2))
boxplot(all_centroids[all_centroids$habitat == habitats[4],]$incl_prob/iters, main = "HSS", ylim = c(0,0.01), ylab = "Inclusion Probability", xlab = paste(habitats[4]))
abline(h = sum(N_Zone[4])/nrow(pts.list[[4]]), col = "red")
boxplot(all_centroids[all_centroids$habitat == habitats[4],]$incl_prob_srs/iters, main = "SRS", ylim = c(0,0.01), ylab = "Inclusion Probability", xlab = paste(habitats[4]))
abline(h = sum(N_Zone[4])/nrow(pts.list[[4]]), col = "red")

# save(all_centroids, file = "./spatial_data/simulated_probs.Rda")