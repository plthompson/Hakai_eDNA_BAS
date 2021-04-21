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
library(cowplot)
library(ggspatial)

source('./code/functions/HSS.r')

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

#create Halton boxes#####
#get bounding box for BC
# MARINE MS, need to add the seed to it for it to work.
bb <- getBB()
attr(bb, "seed") <- getSeed()

box_size <- 25 #chose size of halton box
habitats <- unique(habitat_line_features$habitat)
sampleV <- c(20,20,40,60)
selected_boxes <- list()
HIPBoxes <- list()
for(i in 1:length(habitats)){
  habitat <- habitats[i]
  print(habitat)
  all_boxes <- point2Frame(pts = habitat_line_features %>% filter(habitat == habitats[i]), bb = bb, size = box_size)
  
  coords <- st_coordinates(st_centroid(all_boxes))
  bb.tmp <- st_bbox(bb)
  bbox <- cbind(c(bb.tmp['xmin'], bb.tmp['ymin']), c(bb.tmp['xmax'], bb.tmp['ymax']))
  
  HSS.pts <- getHipSample(X = coords[,1], Y = coords[,2], index = all_boxes$HaltonIndex,
                          N = sampleV[i]*2, bb = bbox,  base = c(2,3), quiet = TRUE,
                          Ps1 = 0:1, Ps2 = 0:2, hipS1 = 0:1, hipS2 = c(0,2,1))	# Changed the order for fun
  n.boxes <- length(table(HSS.pts$HIPIndex))
  
  # Chooses a sample from all boxes. Small wrapper function.
  n <- sampleV[i]
  pts.new <- getSamples(HSS.pts, n = n)
  # Bounding box to clip the HIP boxes to.
  bb.tmp <- st_bbox(research_area)
  bbox2 <- cbind(c(bb.tmp['xmin'], bb.tmp['ymin']), c(bb.tmp['xmax'], bb.tmp['ymax']))
  HIPBoxes[[i]] <- st_as_sf(getHIPBoxes(hip = HSS.pts, bb = bbox2, n = n.boxes, base = c(2,3)))
  HIPBoxes[[i]] <- st_set_crs(HIPBoxes[[i]], value = st_crs(bb))
  
  selected_boxes[[i]] <- all_boxes %>% filter(HaltonIndex %in% pts.new$index) %>% mutate(habitat = habitat)
}

selected_boxes <- rbind(selected_boxes[[1]],selected_boxes[[2]], selected_boxes[[3]], selected_boxes[[4]])

selected_boxes %>% 
  group_by(HaltonIndex) %>% 
  summarise(n = n()) %>% 
  filter(n>1)

selected_boxes_points <- st_centroid(selected_boxes)

area_plot+
  geom_sf(data = habitat_line_features, size = 0.3) +
  geom_sf(data = selected_boxes_points, aes(fill = habitat), color = 1, size = 3, pch = 21)+
  facet_wrap(~habitat)
ggsave("./figures/halton_boxes_facet.pdf", height = 8, width = 8)

area_plot+
  geom_sf(data = habitat_line_features, size = 0.3) +
  geom_sf(data = selected_boxes_points, aes(fill = habitat), color = 1, size = 3, pch = 21)+
  scale_fill_brewer(palette = "Set1", name = "")+
  theme(legend.position = c(1,1), legend.justification = c(1,1))+
  ggspatial::annotation_scale(location = "bl")
ggsave("./figures/halton_boxes.pdf", height = 6, width = 6)

area_plot+
  geom_sf(data = habitat_line_features, size = 0.3) +
  geom_sf(data = selected_boxes, color = 2)+
  scale_fill_brewer(palette = "Set1", name = "")+
  coord_sf(xlim = c(855000,856000), ylim = c(750000, 751000))+
  ggspatial::annotation_scale(location = "bl")
