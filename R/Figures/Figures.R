#making graphs of data###
#this script creates figures of the raw data 

#load libraries
library(h3jsr)
library(raster)
library(sf)
library(patchwork)
library(terra)
library(exactextractr)
library(nlme)
library(tidyverse)

#load species data
plot<-readRDS("traits/species_final_traits.rds")

#load in mean annual temp and temp range 
library(terra)
atr<-terra::rast("climate/wc2.0_bio_10m_07.tif")
mat<-terra::rast("climate/wc2.0_bio_10m_01.tif")

#read in point data for the three states
intersectionsh3<-readRDS("geo_data/points.rds")

#get background maps of the US and a map of the 3 states
states<-readRDS("geo_data/states_polygon.rds")
states2<-readRDS("geo_data/US_polygon.rds")

#using uber h3 to get h3 addresses for resoltution 4
all_res<-point_to_cell(intersectionsh3, res = 4, simple = FALSE)

#quick code here to rename colnames and fix genus names 
colnum<-which(str_detect(colnames(all_res), "resolution"))
colnames(all_res)[colnum]<-"res"
all_res$species<-gsub("\\."," ", all_res$species)

#get exotic % and polymorphic percent for each address
poly_perc<-all_res %>%
  filter(species %in% plot$SPECIES) %>%
  dplyr::select(species,res) %>%
  left_join(.,plot, by = c("species" = "SPECIES")) %>%
  group_by(res) %>%
  filter(!duplicated(species)) %>%
  group_by(res, poly_id) %>%
  summarise(
    S = n()) %>%
  group_by(res)%>%
  mutate(poly_perc = S / sum(S)) %>%
  filter(poly_id == 1) %>%
  ungroup() %>%
  dplyr::select(res, poly_perc)

non_native_perc<-all_res %>%
  filter(species %in% plot$SPECIES) %>%
  dplyr::select(species,res) %>%
  left_join(.,plot, by = c("species" = "SPECIES")) %>%
  group_by(res) %>%
  filter(!duplicated(species)) %>%
  group_by(res, Status) %>%
  summarise(
    S = n()) %>%
  group_by(res)%>%
  mutate(exotic_perc = S / sum(S)) %>%
  filter(Status == "e") %>%
  ungroup() %>%
  dplyr::select(res, exotic_perc)

plotv1<-plot %>%
  mutate(poly_id = ifelse(poly_id == 1, "Polymorphic", "Monomorphic")) %>%
mutate(Status = ifelse(Status == "e", "Nonnative", "Native"))
plotv1$median_size[plotv1$SPECIES == "Dolichoderus.mariae"]<-45250
#plot exotic vs native species
ggplot(data = plotv1 ) +
  geom_point(aes(x = QWD_n, y = log10(median_size),
                 fill = Status), pch = 21,
             size = 4, alpha = 0.8) +
  facet_wrap(~poly_id) +
  geom_smooth(method = "lm",
              aes(x = QWD_n, y = log10(median_size), color = Status))  +
  theme_bw() +
  labs(x = "Queen-worker dimorphism %", y = "Colony Size (Log-scale)") +
  theme(legend.position = "top",
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text( size = 12),
        legend.direction = "horizontal",
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  annotation_logticks(sides = "l")




  
  

###now get average QWD, colony size, and species richness per address
#also filter out cells with less than 5 species
#also remove exotic species
all_res2<-all_res %>%
  filter(species %in% plot$SPECIES) %>%
  dplyr::select(species,res) %>%
  left_join(.,plot, by = c("species" = "SPECIES")) %>%
  filter(!Status == "e") %>%
  group_by(res,poly_id) %>%
  filter(!duplicated(species)) %>%
  summarise(QWD = mean(QWD,na.rm =T),
            QWD_hw = mean(QWD_HW, na.rm =T),
            median_size = median(median_size,na.rm =T),
            var_QWD = var(QWD_HW, na.rm = T),
            S = n()) %>%
  filter(!S<5)

#now get actual hex polygons using all_res2
ant_exo<-cell_to_polygon(all_res2, simple = FALSE) %>%
  st_transform(crs = "EPSG:5070")

#extract mean mat 
se <- cbind(ant_exo, exact_extract(mat, ant_exo, c('mean')))
colnames(se)[which(grepl("exact_extract", colnames(se)))]<-"mean_mat"
#fix poly and mono labels
se<-se %>%
  mutate(poly = ifelse(poly_id == 1, "Polymorphic", "Monomorphic"))

#QWD ~ MAT across morph type
ggplot(data =se) +
  geom_point(aes(x = mean_mat, y = QWD_hw, fill = log10(median_size),
                 size = log10(median_size)), pch = 21, alpha = 0.8) +
  scale_size(guide = 'none') +
  geom_smooth(method = "lm",aes(x = mean_mat, y = QWD_hw)) +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(size = NULL, fill = "Colony size",
       x = "Mean Annual Temp.", y = "Queen-worker Dimorphism") +
  theme(legend.position = "top",
        strip.text =element_text(size = 12))

ggplot(data =se) +
  geom_point(aes(x = mean_mat, y = QWD_hw, fill = as.factor(poly),
                 size = log10(median_size)), pch = 21, alpha = 0.8) +
  scale_size(guide = 'none') +
  geom_smooth(method = "lm",aes(x = mean_mat, y = QWD_hw, color = as.factor(poly))) +
  theme_bw() +
  labs(size = NULL, fill = "Colony size",
       x = "Mean Annual Temp.", y = "Queen-worker Dimorphism") +
  theme(legend.position = "top",
        strip.text =element_text(size = 12))

#average QWD across assemblages without separating out monomorphic and polymorphic 
all_res3<-all_res %>%
  filter(species %in% plot$SPECIES) %>%
  dplyr::select(species,res) %>%
  left_join(.,plot, by = c("species" = "SPECIES")) %>%
  group_by(res) %>%
  filter(!duplicated(species)) %>%
  summarise(QWD = mean(QWD,na.rm =T),
            QWD_hw = mean(QWD_HW, na.rm =T),
            QWD_n = mean(QWD_n, na.rm =T),
            median_size = median(median_size,na.rm =T),
            stdev_QWD = sqrt(var(QWD_HW, na.rm = T)),
            S = n()) %>%
  filter(!S<5)

ant_exo1<-cell_to_polygon(all_res3, simple = FALSE) %>%
  st_transform(crs = "EPSG:5070")
se2 <- cbind(ant_exo1, exact_extract(mat, ant_exo1, c('mean')))
colnames(se)[which(grepl("exact_extract", colnames(se)))]<-"mean_mat"

se2<-left_join(se2, poly_perc, by = c("h3_address" = "res"))
se2<-left_join(se2, non_native_perc, by = c("h3_address" = "res"))

#QWD~ MAT averaged across morph type
both = ggplot(data =se) +
  geom_point(aes(x = mean_mat, y = QWD_hw, fill = log10(median_size),
                 size = log10(median_size)), pch = 21, alpha = 0.8) +
  scale_size(guide = 'none') +
  geom_smooth(method = "lm",aes(x = mean_mat, y = QWD_hw)) +
  scale_fill_viridis_c(option = "B") +
  theme_bw() +
  labs(size = NULL, fill = "Colony size",
       x = "Mean Annual Temp.", y = "Queen-worker Dimorphism",
       title  = "Average of both") +
  theme(legend.position = "top",
        strip.text =element_text(size = 12))


#size ~ QWD with polymorphism percentage
poly = ggplot(data =se2) +
  geom_point(aes(x = QWD_hw, y = log10(median_size), fill = (poly_perc*100),
                 size = (poly_perc*100)), pch = 21, alpha = 0.8) +
  scale_size(guide = 'none') +
  geom_smooth(method = "lm",aes(x = QWD_hw, y = log10(median_size))) +
  scale_fill_viridis_c(option = "C",
                       guide = guide_colorbar(
                         # here comes the code change:
                         frame.colour = "black",
                         ticks = FALSE
                       )) +
  theme_bw() +
  labs(size = NULL, fill = "Assemblage\npolymorphism %",
       x = "Queen-worker Dimorphism %", y = "Colony size (log-scale)") +
  theme(legend.position = c(0.7,0.1),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.title = element_text(face = "bold", size = 12),
        legend.direction = "horizontal") +
  annotation_logticks(sides = "l")
poly


#QWD variation vs QWD average
ggplot(data = se2) +
  geom_point(aes(x = stdev_QWD, y = QWD_hw))

#map of QWD
poly_map<-ggplot() +
  geom_sf(data = states2) +
  geom_sf(data = se2, aes(fill = QWD_n, color = QWD_n)) +
  scale_fill_viridis_c(option = "D",
                       guide = guide_colorbar(
                         # here comes the code change:
                         frame.colour = "black",
                         ticks = FALSE
                       )) +
  scale_color_viridis_c(option = "D") +
  theme_bw() +
  coord_sf(xlim = c(605755.4,
                    1632548.473),
           ylim = c(
             1406449.988, 250506.46)) +
  theme(panel.background = element_rect(fill = "skyblue"),
        legend.position = c(0.3,0.1),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.title = element_text(face = "bold", size = 12),
        legend.direction = "horizontal") +
  labs(fill = "Queen-worker\ndimorphism %", color = "Queen-worker\ndimorphism %")
poly_map

poly_map + poly


#variance ~ mean
ggplot(data =se2) +
  geom_point(aes(x = stdev_QWD, y = QWD_hw, fill = (poly_perc),
                 size = (poly_perc)), pch = 21, alpha = 0.8) +
  scale_size(guide = 'none') +
  geom_smooth(method = "lm",aes(x = stdev_QWD, y = QWD_hw))


#also now get species-averages
sp_avg<-cbind(terra::extract(mat,intersectionsh3), intersectionsh3)
sp_avg$species<-gsub("\\."," ", sp_avg$species)
sp_avg$wc2.0_bio_10m_07
sp_avg1<-sp_avg %>%
  group_by(species) %>%
  summarise(environ = sqrt(var(wc2.0_bio_10m_01)),
            n = n()) %>%
  filter(n >20) %>%
  left_join(., plot, by = c( "species" = "SPECIES")) %>%
  drop_na(poly_id)

ggplot(sp_avg1) +
  geom_point(aes(x = environ,
                 y = (QWD_n), color = as.factor(Status))) +
  geom_smooth( method = "lm",aes(x = environ, y = QWD_HW,
                                 color = as.factor(Status)) ) +
  facet_wrap(~poly_id)


#compare species richness of raster to raw data

