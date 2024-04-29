#urban data extraction
library(h3jsr)
library(sf)
library(terra)
library(landscapemetrics)
library(exactextractr)
library(tidyverse)
#load CONUS LC data
#Now load in your data from the PGLS resampling
resampled1<-readRDS("pgls_results/FBD_crown_resampled_nonnative_concensus_res4.rds")

#read in the original data that you used to generate pgls results
pgls_data<-readRDS("pgls_results/PGLS_datasets_resample_nonnative/pgls_data_res4.rds")

#read in native percentage
native_perc<-readRDS("pgls_results/exotic_percentage_forpgls_res4.rds")

resampled1<-resampled1 %>%
  filter(!assemblage == "null") %>%
  filter(homogeneity > 0.05) %>%
  filter(normality >0.05) %>%
  filter(!r2 == "error") %>%
  filter(r2 <1 & r2 > 0)


resampled1<-left_join(resampled1, native_perc, by = c("ID", "assemblage"))


#gives assemblage based slopes, intercepts, pvals etc.
resampled<-resampled1  %>%
  group_by(res) %>%
  mutate_all(~ as.numeric(as.character(.))) %>%
  summarise(slope = mean(slope, na.rm  =T),
            intercept = mean(intercept, na.rm = T),
            slope_err = (mean(slope_error, na.rm =T)),
            avg_r2 = mean(r2),
            pval = median(p)) 



#get hexagons for the data
hexes = readRDS("geo_data/non_native_hexes_figures_res4.rds")
pgls_data1<-pgls_data %>%
  dplyr::select(res, ID)

#Give every passing model the res h3 code.
resampled2<-left_join(resampled1, pgls_data1, by = "res")

#give assemblage average their res codes
resampled1<-left_join(resampled, pgls_data, by = "res")

#join with lc data from landscapemetrics
urbanhmetrics1<-readRDS("geo_data/lanscapemetrics_nonnative_lc_res4.rds")
wtf<-left_join(urbanhmetrics1, resampled1, by = c("res" = "res"))

wtf2<-wtf %>% filter(pval<0.05)

#remove unnneeded classes 
wtf2<-wtf2 %>% filter(class %in% c(1,2,3,4)) %>%
  mutate(class = case_when(
    class == 1 ~ "Developed",
    class == 2 ~ "Cropland",
    class == 3 ~ "Grass_Shrub",
    class == 4 ~ "Tree_Cover"
  ))

diff_metrics<-wtf2 %>%
  pivot_wider(id_cols = c(res,class), names_from = metric, values_from = value)


saveRDS(diff_metrics, "native_nonnative_results/landscape_metrics_nonnative_res4.rds")
saveRDS(resampled2, "native_nonnative_results/nonnative_landscape_results_res4.rds")


urbanhmetrics_area<-readRDS("geo_data/lanscapemetrics_nonnative_lc_area_res4.rds")
#convert from hectares to km2
urbanhmetrics_area$value<-urbanhmetrics_area$value/100
urbanhmetrics_area1<-urbanhmetrics_area %>%
  filter(value >200) %>%
  dplyr::select(res, value) %>%
  rename(area = value)
saveRDS(urbanhmetrics_area1, "native_nonnative_results/area_corrected_nonnative_res4.rds")


metrics_GI<-readRDS("geo_data/lanscapemetrics_nonnative_lc_GI_res4.rds")

metrics_GI1<-metrics_GI %>%
  dplyr::select(res, value) %>%
  rename(Conditional_entropy = value)


#remove areas where we have low area coverage
wtf3<-wtf2 %>%
  filter(res %in% urbanhmetrics_area1$res) %>%
  dplyr::select(res, class, metric, value,slope, avg_r2, pval,intercept) %>%
  pivot_wider(id_cols = c(res, class, slope,avg_r2,intercept), names_from = metric, values_from = value)

#get coordinates of ids
wtf3 <-left_join(wtf3,metrics_GI1, by = c("res" = "res"))

coordinates<-cell_to_polygon(pgls_data1$res) %>%
  st_as_sf() %>%
  st_transform(crs = 4326) %>%
  st_centroid() %>%
  st_coordinates() %>%
  as.data.frame()

res<-data.frame(res = pgls_data1$res)

xy<-cbind(res, coordinates)

landscape_metrics_class<-left_join(wtf3, xy, by = "res") %>%
  pivot_wider(id_cols = c(res,X,Y), 
              names_from = class, 
              values_from = c(Total_class_area, Clumpiness_index, Contiguity_index))

ID_row<-unique(wtf3[, c("res","slope", "avg_r2","intercept")])
model_data_nonnative<-left_join(ID_row, landscape_metrics_class, by = "res")
saveRDS(model_data_nonnative, "geo_data/spatial_glm_data/nonnative/model_all_res4.rds")






