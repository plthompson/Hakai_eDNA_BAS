# code by Patrick Thompson
# April 2021
# patrick.thompson@dfo-mpo.gc.ca

#load packages#####


#library(usethis)
#library(devtools)
#devtools::install_github("paul-vdb/DFO-master-sample")

library(BASMasterSample)
library(sp)
library(rgdal)
library(sf)
library(tidyverse)
library(cowplot)
library(RColorBrewer)

#ben's directory
setwd("C:/Users/Ben/OneDrive - McGill University/Documents/eDNA Thesis/repos/Hakai_eDNA_BAS")

#change plot theme####
source("./code/functions/plot_theme.R")
theme_set(plt_theme)

#load research area####
research_area <- readOGR("./spatial_data/calvert_research_area.gpkg")
plot(research_area)
research_area <- st_as_sf(research_area, wtk = geometry)
area_plot <- ggplot()+
  geom_sf(data = research_area, fill = NA)

area_plot

end
#polygon based BAS####
# this is equivalent to the random method that Ben used, but with two advantages
# 1 - it should have better spatial balance - I need to check this because the polygons don't represent a continuous area
# 2 - the samples have an order and this could be used to ask how much spatial coverage is needed before biodiversity saturates

#read in habitat polygons
habitat_polygon_features <- readOGR("./spatial_data/habitat_polygon_features_discretized.gpkg")
habitat_polygon_features <- st_as_sf(habitat_polygon_features)
habitat_polygon_features
#this shapefile has legacy buffers around sites and "no-go" zones removed


#read in legacy sites
nearshore_legacy_sites <- readOGR("./spatial_data/hakai_legacy_sites.gpkg")

nearshore_legacy_sites <- st_as_sf(nearshore_legacy_sites)

#stratified BAS sample
#choose how many sites in each polygon in the dataset - I matched these to the numbers in Ben's proposal
# legacy sites to subtract: bull_kelp 5, giant_kelp 7, high_rugosity 0, low_rugosity 0, seagrass 16, unclassified 14

N_Zone <- c("high_rugosity" = 40,
            "low_rugosity" = 20,
            "bull_kelp" = 15,
            "giant_kelp" = 13,
            "seagrass" = 24,
            "unclassified" = 46)

areaBAS <- masterSample(shp = habitat_polygon_features, N = N_Zone, stratum = "habitat")

areaBAS$habitat <- c(rep(names(N_Zone)[1], N_Zone[1]),
                   rep(names(N_Zone)[2], N_Zone[2]),
                   rep(names(N_Zone)[3], N_Zone[3]),
                   rep(names(N_Zone)[4], N_Zone[4]),
                   rep(names(N_Zone)[5], N_Zone[5]),
                   rep(names(N_Zone)[6], N_Zone[6]))


table(areaBAS$SiteID)
#how many unique sample locations do we have? (because some overlap across the different habitat layers)
areaBAS %>%
  select(-habitat) %>%
  unique() %>%
  nrow()

#show sites that are included in multiple layers
areaBAS <- st_join(areaBAS, areaBAS[areaBAS %>%
                                      select(-layer) %>%
                                      duplicated(),] %>%
                     select(SiteID) %>%
                     mutate(multiple = TRUE)) %>%
  select(-SiteID.y, SiteID = SiteID.x)

areaBAS$multiple[is.na(areaBAS$multiple)] <- FALSE

coordinates <- data.frame(st_coordinates(areaBAS))
coordinates$SiteID <- areaBAS$SiteID
coordinates$Layer <- areaBAS$layer

ggplot(coordinates, aes(x = SiteID, y = Y))+
  geom_point()+
  facet_wrap(~Layer, scales = "free")

ggplot(coordinates, aes(x = SiteID, y = X))+
  geom_point()+
  facet_wrap(~Layer, scales = "free")

ggplot(coordinates, aes(x = X, y = Y, color = SiteID))+
  geom_path()+
  facet_wrap(~Layer, scales = "free")+
  scale_color_viridis_c()

ggplot()+
  geom_sf(data = areaBAS, aes(fill = layer), size = 2, pch = 21)+
  scale_fill_brewer(palette = "Set1", name = "")

ggplot()+
  geom_sf(data = research_area, fill = NA)+
  geom_sf(data = habitat_polygon_features, size = 0.1)+
  geom_sf(data = areaBAS, aes(fill = habitat), size = 2, pch = 21)+
  geom_sf(data = nearshore_legacy_sites, aes(fill = habitat), size = 2, pch = 21)+
  scale_fill_brewer(palette = "Set1", name = "")
#ggsave("./figures/BAS_polygon_method.pdf", height = 6, width = 6)

ggplot()+
  geom_sf(data = research_area, fill = NA)+
  geom_sf(data = nearshore_legacy_sites, aes(fill = Habitat), size = 2, pch = 21)+
  scale_fill_brewer(palette = "Set1", name = "")
ggsave("./figures/legacy sites.pdf", height = 6, width = 6)

writeOGR(areaBAS, "areaBAS", driver="ESRI Shapefile")

?writeOGR

orgDrivers()

#halton box method#####
# note - bounding box should be all of BC to make boxes align with greater BC BAS
# however this is slow and so for now we will use the study region as the bounding box
#make bounding box
bb <- buildMS(research_area, d = 9, FALSE)

#first check that halton boxes of different size line up
halton_boxes_0.02 <- point2Frame(pts = habitat_line_features, bb = bb, size = 1000) %>%
  arrange(HaltonIndex) %>%
  top_n(wt = HaltonIndex, n = -4) #select the number of boxes. The - sign indicates to take the lowest values

halton_boxes_0.05 <- point2Frame(pts = habitat_line_features, bb = bb, size = 2000) %>%
  arrange(HaltonIndex) %>%
  top_n(wt = HaltonIndex, n = -6)

# smaller boxes always fall within larger boxes
# larger boxes only contain smaller boxes if habitat is present within the subset of the box encompased by the smaller box
area_plot_color+
  geom_sf(data = halton_boxes_0.05, fill = "red")+
  geom_sf(data = halton_boxes_0.02, fill = "black")


#select halton boxes for each habitat type
boxes <- 50 #choose number of sites to include
box_size <- 0.02 #chose size of halton box

halton_boxes_seagrass <- point2Frame(pts = habitat_polygon_features_legacyremoved %>% filter(hab == "seagrass"), bb = bb, size = 1000) %>%
  arrange(HaltonIndex) %>%
  top_n(wt = HaltonIndex, n = -40)

halton_boxes_bull_kelp <- point2Frame(pts = habitat_polygon_features_legacyremoved %>% filter(hab == "bull_kelp"), bb = bb, size = 1000) %>%
  arrange(HaltonIndex) %>%
  top_n(wt = HaltonIndex, n = -20)

halton_boxes_giant_kelp <- point2Frame(pts = habitat_polygon_features_legacyremoved %>% filter(hab == "giant_kelp"), bb = bb, size = 1000) %>%
  arrange(HaltonIndex) %>%
  top_n(wt = HaltonIndex, n = -20)

halton_boxes_unclassified <- point2Frame(pts = habitat_polygon_features_legacyremoved %>% filter(hab == "unclassified"), bb = bb, size = 1000) %>%
  arrange(HaltonIndex) %>%
  top_n(wt = HaltonIndex, n = -60)

halton_boxes_highrugosity <- point2Frame(pts = habitat_polygon_features_legacyremoved %>% filter(hab == "high_rugosity"), bb = bb, size = 1000) %>%
  arrange(HaltonIndex) %>%
  top_n(wt = HaltonIndex, n = -40)

halton_boxes_lowrugosity <- point2Frame(pts = habitat_polygon_features_legacyremoved %>% filter(hab == "low_rugosity"), bb = bb, size = 1000) %>%
  arrange(HaltonIndex) %>%
  top_n(wt = HaltonIndex, n = -20)




sample.df <- full_join(data.frame(HaltonIndex = halton_boxes_seagrass$HaltonIndex,Seagrass = TRUE),
                       data.frame(HaltonIndex = halton_boxes_bull_kelp$HaltonIndex, BullKelp = TRUE)) %>%
  full_join(data.frame(HaltonIndex = halton_boxes_giant_kelp$HaltonIndex, GiantKelp = TRUE)) %>%
  full_join(data.frame(HaltonIndex = halton_boxes_unclassified$HaltonIndex, Unclassified = TRUE)) %>%
  full_join(data.frame(HaltonIndex = halton_boxes_highrugosity$HaltonIndex, Unclassified = TRUE)) %>%
  full_join(data.frame(HaltonIndex = halton_boxes_lowrugosity$HaltonIndex, Unclassified = TRUE)) %>%
  arrange(HaltonIndex) %>%
  group_by(HaltonIndex) %>%
  mutate(Habitats = sum(Seagrass, BullKelp, GiantKelp, Unclassified, na.rm = TRUE))



ggplot(sample.df, aes(x = Habitats))+
  geom_histogram(bins = 4)

ggplot(sample.df, aes(x = HaltonIndex, y = Habitats))+
  geom_point()

colV <- brewer.pal(4, name = "Dark2")

plot_grid(
  area_plot_color+
    geom_sf(data = halton_boxes_bull_kelp, fill = NA, size = 0.3, color = 1)+
    ggtitle(label = "Bull kelp"),

  area_plot_color+
    geom_sf(data = halton_boxes_giant_kelp, fill = NA, size = 0.3, color = 1)+
    ggtitle(label = "Giant kelp"),

  area_plot_color+
    geom_sf(data = halton_boxes_seagrass, fill = NA, size = 0.3, color = 1)+
    ggtitle(label = "Seagrass"),

  area_plot_color+
    geom_sf(data = halton_boxes_unclassified, fill = NA, size = 0.3, color = 1)+
    ggtitle(label = "Unclassified"),
  nrow = 2)


area_plot_color+
  geom_sf(data = halton_boxes_highrugosity, fill = NA, size = 0.3, color = 1)+
  ggtitle(label = "high rugosity")

area_plot_color+
  geom_sf(data = halton_boxes_lowrugosity, fill = NA, size = 0.3, color = 1)+
  ggtitle(label = "high rugosity")

#ggsave("./figures/halton_boxes.pdf", height = 8, width = 6)

#plot all together
area_plot+
  geom_sf(data = halton_boxes_unclassified, fill = colV[4], alpha = 1)+
  geom_sf(data = halton_boxes_giant_kelp, fill = colV[2], alpha = 1)+
  geom_sf(data = halton_boxes_bull_kelp, fill = colV[1], alpha = 1)+
  geom_sf(data = halton_boxes_seagrass, fill = colV[3], alpha = 1)+
  geom_sf(data = halton_boxes_highrugosity, fill = colV[1], alpha = 1)+
  geom_sf(data = halton_boxes_lowrugosity, fill = colV[3], alpha = 1)+
  geom_sf(data = nearshore_legacy_sites, aes(color = 'blue'), alpha = 0.3)

