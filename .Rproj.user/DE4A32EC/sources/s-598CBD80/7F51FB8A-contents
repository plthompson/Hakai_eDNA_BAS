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

#read in habitat lines
habitat_line_features <- readOGR("./spatial_data/habitat_line_features.gpkg")
habitat_line_features <- st_as_sf(habitat_line_features)

#read in habitat polygons
habitat_polygons <- readOGR("./spatial_data/habitat_line_polyRug_legbuff.gpkg")
habitat_polygons <- st_as_sf(habitat_polygons)
plot(habitat_polygons)

#function to produce HIP BAS for a habitat
#habitats is the habitat that you are selecting sites for
#priority habitats are habitats that would be preferentially sampled in a given halton box
#this excludes any halton box with those habitats from consideration so that only boxes without priority habitats are selected
habitat_HIP <- function(habitats, samples, priority_habitats = NA){
  #create Halton boxes#####
  #get bounding box for BC
  # MARINE MS, need to add the seed to it for it to work.
  bb <- getBB()
  attr(bb, "seed") <- getSeed()
  
  all_boxes <- point2Frame(pts = habitat_polygons %>% filter(habitat %in% habitats), bb = bb, size = 100)
  if(!is.na(priority_habitats)){
    all_boxes_priority <- point2Frame(pts = habitat_polygons %>% filter(habitat %in% priority_habitats), bb = bb, size = 100)
    all_boxes <- all_boxes %>% filter(!(HaltonIndex %in% all_boxes_priority$HaltonIndex))
  }
  
  coords <- st_coordinates(st_centroid(all_boxes))
  bb.tmp <- st_bbox(bb)
  bbox <- cbind(c(bb.tmp['xmin'], bb.tmp['ymin']), c(bb.tmp['xmax'], bb.tmp['ymax']))
  
  HSS.pts <- getHipSample(X = coords[,1], Y = coords[,2], index = all_boxes$HaltonIndex,
                          N = samples*2, bb = bbox,  base = c(2,3), quiet = TRUE,
                          Ps1 = 0:1, Ps2 = 0:2, hipS1 = 0:1, hipS2 = c(0,2,1))	# Changed the order for fun
  n.boxes <- length(table(HSS.pts$HIPIndex))
  
  # Chooses a sample from all boxes. Small wrapper function.
  pts.new <- getSamples(HSS.pts, n = samples)
  # Bounding box to clip the HIP boxes to.
  bb.tmp <- st_bbox(research_area)
  bbox2 <- cbind(c(bb.tmp['xmin'], bb.tmp['ymin']), c(bb.tmp['xmax'], bb.tmp['ymax']))
  HIPBoxes <- st_as_sf(getHIPBoxes(hip = HSS.pts, bb = bbox2, n = n.boxes, base = c(2,3)))
  HIPBoxes <- st_set_crs(HIPBoxes, value = st_crs(bb))
  
  plot(HIPBoxes)
  
  selected_boxes <- all_boxes %>% filter(HaltonIndex %in% pts.new$index)
  
  selected_boxes_points <- st_centroid(selected_boxes)
  
  area_plot+
    geom_sf(data = habitat_line_features, size = 0.3) +
    geom_sf(data = selected_boxes_points, color = 2, size = 3)
  
  selected_boxes <- st_intersection(habitat_polygons %>%  filter(habitat %in% habitats), selected_boxes)
  return(selected_boxes)
}


bb <- getBB()
attr(bb, "seed") <- getSeed()
  


#get samples per habitat####
unclassified_sites <- habitat_HIP(habitats = "unclassified", samples = 40, priority_habitats = c("giant_kelp", "bull_kelp", "seagrass"))
deeper_sites <- masterSample(shp = habitat_polygons %>% filter(habitat %in% c("low_rugosity", "high_rugosity")), N = c("low_rugosity" = 20, "high_rugosity" = 40), stratum = "habitat", bb = bb)
deeper_sites$habitat <- c(rep("low_rugosity", 20), rep("high_rugosity", 40))
#these next two are if you want to use halton boxes for the low and high rugosity sites
#low_rugosity_sites <- habitat_HIP(habitats = "low_rugosity", samples = 20, priority_habitats = c("giant_kelp", "bull_kelp", "seagrass", "unclassified", "high_rugosity"))
#high_rugosity_sites <- habitat_HIP(habitats = "high_rugosity", samples = 40, priority_habitats = c("giant_kelp", "bull_kelp", "seagrass", "unclassified"))
habitat_sites <- habitat_HIP(habitats = c("giant_kelp", "bull_kelp", "seagrass"), samples = 68)


all_sites <- bind_rows(unclassified_sites, habitat_sites) %>% #add low and high rugosity sites here if you want to use the halton box method for them
  group_by(HaltonIndex) %>% 
  mutate(count = n()) %>% 
  arrange(desc(count), HaltonIndex)

selected_boxes_points <- st_centroid(all_sites)

box_count <- all_sites %>% 
  as.data.frame() %>% 
  dplyr::select(HaltonIndex, count) %>% 
  unique()

table(box_count$count)
table(all_sites$habitat) #note this includes overlapping sites for classified habitat so this number will be lower in reality

area_plot+
  geom_sf(data = habitat_line_features, size = 0.3) +
  geom_sf(data = selected_boxes_points, aes(fill = factor(count)), size = 3, pch = 21)+
  geom_sf(data = deeper_sites, size = 3, pch = 21, fill = "#e41a1c")+
  scale_fill_brewer(palette = "Set1", name = "habitats\nin box")+
  facet_wrap(~habitat)
ggsave("./figures/BAS_general_2.pdf", height = 8, width = 12)

