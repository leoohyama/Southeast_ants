#calculate landscape metrics for all h3 cells from res 5 to 3

library(sf)
library(terra)
library(h3jsr)
library(tidyverse)
library(landscapemetrics)
#load CONUS LC data
states_raster<-terra::rast("geo_data/state_raster_lc.tif")

plot(states_raster)
#get a polygon of the states
polygon_states<-readRDS("geo_data/states_polygon.rds")
#read in polygon of states based on lc data

polygon_states_buffer<-st_buffer(polygon_states, dist = 60000)

polygon_states_reproj =polygon_states %>%
  st_transform(crs = st_crs(states_raster))

#get addresses for all three states using uber h3 and generate polygons
res4_address<-polygon_to_cells(polygon_states_buffer, res =4)
res5_address<-polygon_to_cells(polygon_states_buffer, res =5)
res3_address<-polygon_to_cells(polygon_states_buffer, res =3)
res6_address<-polygon_to_cells(polygon_states_buffer, res =6)

#get hexagons for the data
hexes_res4 = cell_to_polygon(res4_address, simple = FALSE) %>%
  st_transform(crs = st_crs(states_raster))
hexes_res5 = cell_to_polygon(res5_address, simple = FALSE)%>%
  st_transform(crs = st_crs(states_raster))
hexes_res3 = cell_to_polygon(res3_address, simple = FALSE)%>%
  st_transform(crs = st_crs(states_raster))
hexes_res6 = cell_to_polygon(res6_address, simple = FALSE)%>%
  st_transform(crs = st_crs(states_raster))

#crop all res 4 to the landcover data
polygon_states_lc<-readRDS("geo_data/LC_polygon_allstates.rds")
intersected_res4<-st_intersection(hexes_res4, polygon_states_lc)
intersected_res3<-st_intersection(hexes_res3, polygon_states_lc)
intersected_res5<-st_intersection(hexes_res5, polygon_states_lc)
intersected_res6<-st_intersection(hexes_res6, polygon_states_lc)



directory_lc="geo_data/hex_rasters_lc_res4/"

cropped_hexes<-lapply(1:nrow(intersected_res4), function(x){
  print(x)
  cropped<-terra::crop(states_raster,vect(intersected_res4$geometry[x]), mask = T)
  filename<-paste(intersected_res4$h3_address[x], "lc.tif", sep = "_")
  path_lc=paste(directory_lc,filename,sep = "/")
  writeRaster(cropped, path_lc, overwrite =T)
  return(cropped)
})


directory_lc="geo_data/hex_rasters_lc_res3/"
cropped_hexes_res3<-lapply(1:nrow(intersected_res3), function(x){
  cropped<-terra::crop(states_raster,vect(intersected_res3$geometry[x]), mask = T)
  filename<-paste(intersected_res3$h3_address[x], "lc.tif", sep = "_")
  path_lc=paste(directory_lc,filename,sep = "/")
  writeRaster(cropped, path_lc, overwrite =T)
  return(cropped)
})

directory_lc="geo_data/hex_rasters_lc_res5/"
cropped_hexes_res5<-lapply(1:nrow(intersected_res5), function(x){
  cropped<-terra::crop(states_raster,vect(intersected_res5$geometry[x]), mask = T)
  filename<-paste(intersected_res5$h3_address[x], "lc.tif", sep = "_")
  path_lc=paste(directory_lc,filename,sep = "/")
  writeRaster(cropped, path_lc, overwrite =T)
  return(cropped)
})

directory_lc="geo_data/hex_rasters_lc_res6/"
cropped_hexes_res6<-lapply(1:nrow(intersected_res6), function(x){
  cropped<-terra::crop(states_raster,vect(intersected_res6$geometry[x]), mask = T)
  filename<-paste(intersected_res6$h3_address[x], "lc.tif", sep = "_")
  path_lc=paste(directory_lc,filename,sep = "/")
  writeRaster(cropped, path_lc, overwrite =T)
  return(cropped)
})

