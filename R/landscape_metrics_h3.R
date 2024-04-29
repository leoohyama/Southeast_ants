#calculate landscape metrics for all h3 cells from res 5 to 3
library(tidyverse) 
library(stringr)
library(sf)
library(terra)
library(h3jsr)
library(tidyverse)
library(landscapemetrics)


#get a polygon of the states
polygon_states<-readRDS("geo_data/states_polygon.rds")
plot(polygon_states)

#get addresses for all three states using uber h3 and generate polygons
res4_address<-polygon_to_cells(polygon_states, res =4)
res5_address<-polygon_to_cells(polygon_states, res =5)
res3_address<-polygon_to_cells(polygon_states, res =3)
res6_address<-polygon_to_cells(polygon_states, res =6)

#get hexagons for the data
hexes_res4 = cell_to_polygon(res4_address, simple = FALSE) %>%
  st_transform(crs = st_crs(states_raster))
hexes_res5 = cell_to_polygon(res5_address, simple = FALSE)%>%
  st_transform(crs = st_crs(states_raster))
hexes_res3 = cell_to_polygon(res3_address, simple = FALSE)%>%
  st_transform(crs = st_crs(states_raster))
hexes_res6 = cell_to_polygon(res6_address, simple = FALSE)%>%
  st_transform(crs = st_crs(states_raster))


directory_lc_res4="geo_data/hex_rasters_lc_res4/"
directory_lc_res3="geo_data/hex_rasters_lc_res3/"
directory_lc_res5="geo_data/hex_rasters_lc_res5/"
directory_lc_res6="geo_data/hex_rasters_lc_res6/"

res4_files<-list.files(directory_lc_res4)
res5_files<-list.files(directory_lc_res5)
res6_files<-list.files(directory_lc_res6)
res3_files<-list.files(directory_lc_res3)

#get metrics for every hexbin

#res4
hex_area_res4<-lapply(seq_along(res4_files), function(x){
  raster_read<-paste(directory_lc_res4, res4_files[[x]], sep = "/")
  raster_calc<-terra::rast(raster_read)
  calcs<-tryCatch(
    {calculate_lsm(raster_calc, what = c("lsm_l_ta"))},
    error = function(err) {return("error")}
  )
  
  rm(raster_calc)
  
  if(length(calcs)==1){ return(calcs)}else{
    calcs$res =str_extract(res4_files[[x]], "[^_]+")
    return(calcs)
  }
})

conditonal_entropy_res4<-lapply(seq_along(res4_files), function(x){
  raster_read<-paste(directory_lc_res4, res4_files[[x]], sep = "/")
  raster_calc<-terra::rast(raster_read)
  calcs<-tryCatch(
    {calculate_lsm(raster_calc, what = c("lsm_l_condent"))},
    error = function(err) {return("error")}
  )
  
  rm(raster_calc)
  
  if(length(calcs)==1){ return(calcs)}else{
    calcs$res =str_extract(res4_files[[x]], "[^_]+")
    return(calcs)
  }
})

hex_metrics_res4<-lapply(seq_along(res4_files), function(x){
  raster_read<-paste(directory_lc_res4, res4_files[[x]], sep = "/")
  raster_calc<-terra::rast(raster_read)
  plot(raster_calc)
  calcs<-tryCatch(
    {calculate_lsm(raster_calc, what = c("lsm_c_clumpy","lsm_c_ca",
                                         "lsm_c_contig_mn"))},
    
    error = function(err) {return("error")}
  )

  rm(raster_calc)
  
  if(length(calcs)==1){ return(calcs)}else{
    
    calcs$res = str_extract(res4_files[[x]], "[^_]+")
    return(calcs)
  }
  
})


remove<-which(hex_metrics_res4 == "error")
if(length(remove) == 0){
  conditonal_entropy_res4<-conditonal_entropy_res4
  hex_area_res4<-hex_area_res4
  hex_metrics_res4<-hex_metrics_res4
}else{
  conditonal_entropy_res4<-conditonal_entropy_res4[-remove]
  hex_area_res4<-hex_area_res4[-remove]
  hex_metrics_res4<-hex_metrics_res4[-remove]
}



ce_res4<-do.call(rbind, conditonal_entropy_res4) %>%
  dplyr::select(res, value) %>%
  rename(condent = value)
area_res4<-do.call(rbind, hex_area_res4) %>%
  dplyr::select(res, value)
metrics_res4<-do.call(rbind, hex_metrics_res4) %>%
  dplyr::select(res,class,metric, value)


urbanhmetrics1<-left_join(metrics_res4, ce_res4, by = "res") %>%
  mutate(metric = case_when(
    metric == "clumpy" ~ "Clumpiness_index",
    metric == "ca" ~ "Total_class_area",
    metric == "contig_mn" ~ "Contiguity_index"
  ))



saveRDS(urbanhmetrics1, "geo_data/lc_res4.rds")
saveRDS(area_res4, "geo_data/lc_area_res4.rds")
saveRDS(ce_res4, "geo_data/lc_condent_res4.rds")




#res6
hex_area_res6<-lapply(seq_along(res6_files), function(x){
  raster_read<-paste(directory_lc_res6, res6_files[[x]], sep = "/")
  raster_calc<-terra::rast(raster_read)
  calcs<-tryCatch(
    {calculate_lsm(raster_calc, what = c("lsm_l_ta"))},
    error = function(err) {return("error")}
  )
  
  rm(raster_calc)
  
  if(length(calcs)==1){ return(calcs)}else{
    calcs$res =str_extract(res6_files[[x]], "[^_]+")
    return(calcs)
  }
})

conditonal_entropy_res6<-lapply(seq_along(res6_files), function(x){
  raster_read<-paste(directory_lc_res6, res6_files[[x]], sep = "/")
  raster_calc<-terra::rast(raster_read)
  calcs<-tryCatch(
    {calculate_lsm(raster_calc, what = c("lsm_l_condent"))},
    error = function(err) {return("error")}
  )
  
  rm(raster_calc)
  
  if(length(calcs)==1){ return(calcs)}else{
    calcs$res =str_extract(res6_files[[x]], "[^_]+")
    return(calcs)
  }
})

hex_metrics_res6<-lapply(seq_along(res6_files), function(x){
  raster_read<-paste(directory_lc_res6, res6_files[[x]], sep = "/")
  raster_calc<-terra::rast(raster_read)
  plot(raster_calc)
  calcs<-tryCatch(
    {calculate_lsm(raster_calc, what = c("lsm_c_clumpy","lsm_c_ca",
                                         "lsm_c_contig_mn"))},
    
    error = function(err) {return("error")}
  )
  
  rm(raster_calc)
  
  if(length(calcs)==1){ return(calcs)}else{
    
    calcs$res = str_extract(res6_files[[x]], "[^_]+")
    return(calcs)
  }
  
})
remove<-which(hex_metrics_res6 == "error")
if(length(remove) == 0){
  conditonal_entropy_res6<-conditonal_entropy_res6
  hex_area_res6<-hex_area_res6
  hex_metrics_res6<-hex_metrics_res6
}else{
  conditonal_entropy_res6<-conditonal_entropy_res6[-remove]
  hex_area_res6<-hex_area_res6[-remove]
  hex_metrics_res6<-hex_metrics_res6[-remove]
}



ce_res6<-do.call(rbind, conditonal_entropy_res6) %>%
  dplyr::select(res, value) %>%
  rename(condent = value)
area_res6<-do.call(rbind, hex_area_res6) %>%
  dplyr::select(res, value)
metrics_res6<-do.call(rbind, hex_metrics_res6) %>%
  dplyr::select(res,class,metric, value)


urbanhmetrics1<-left_join(metrics_res6, ce_res6, by = "res") %>%
  mutate(metric = case_when(
    metric == "clumpy" ~ "Clumpiness_index",
    metric == "ca" ~ "Total_class_area",
    metric == "contig_mn" ~ "Contiguity_index"
  ))



saveRDS(urbanhmetrics1, "geo_data/lc_res6.rds")
saveRDS(area_res6, "geo_data/lc_area_res6.rds")
saveRDS(ce_res6, "geo_data/lc_condent_res6.rds")



#res5
hex_area_res5<-lapply(seq_along(res5_files), function(x){
  raster_read<-paste(directory_lc_res5, res5_files[[x]], sep = "/")
  raster_calc<-terra::rast(raster_read)
  calcs<-tryCatch(
    {calculate_lsm(raster_calc, what = c("lsm_l_ta"))},
    error = function(err) {return("error")}
  )
  
  rm(raster_calc)
  
  if(length(calcs)==1){ return(calcs)}else{
    calcs$res =str_extract(res5_files[[x]], "[^_]+")
    return(calcs)
  }
})

conditonal_entropy_res5<-lapply(seq_along(res5_files), function(x){
  raster_read<-paste(directory_lc_res5, res5_files[[x]], sep = "/")
  raster_calc<-terra::rast(raster_read)
  calcs<-tryCatch(
    {calculate_lsm(raster_calc, what = c("lsm_l_condent"))},
    error = function(err) {return("error")}
  )
  
  rm(raster_calc)
  
  if(length(calcs)==1){ return(calcs)}else{
    calcs$res =str_extract(res5_files[[x]], "[^_]+")
    return(calcs)
  }
})

hex_metrics_res5<-lapply(seq_along(res5_files), function(x){
  raster_read<-paste(directory_lc_res5, res5_files[[x]], sep = "/")
  raster_calc<-terra::rast(raster_read)
  plot(raster_calc)
  calcs<-tryCatch(
    {calculate_lsm(raster_calc, what = c("lsm_c_clumpy","lsm_c_ca",
                                         "lsm_c_contig_mn"))},
    
    error = function(err) {return("error")}
  )
  
  rm(raster_calc)
  
  if(length(calcs)==1){ return(calcs)}else{
    
    calcs$res = str_extract(res5_files[[x]], "[^_]+")
    return(calcs)
  }
  
})

remove<-which(hex_metrics_res5 == "error")
if(length(remove) == 0){
  conditonal_entropy_res5<-conditonal_entropy_res5
  hex_area_res5<-hex_area_res5
  hex_metrics_res5<-hex_metrics_res5
}else{
  conditonal_entropy_res5<-conditonal_entropy_res5[-remove]
  hex_area_res5<-hex_area_res5[-remove]
  hex_metrics_res5<-hex_metrics_res5[-remove]
}



ce_res5<-do.call(rbind, conditonal_entropy_res5) %>%
  dplyr::select(res, value) %>%
  rename(condent = value)
area_res5<-do.call(rbind, hex_area_res5) %>%
  dplyr::select(res, value)
metrics_res5<-do.call(rbind, hex_metrics_res5) %>%
  dplyr::select(res,class,metric, value)


urbanhmetrics1<-left_join(metrics_res5, ce_res5, by = "res") %>%
  mutate(metric = case_when(
    metric == "clumpy" ~ "Clumpiness_index",
    metric == "ca" ~ "Total_class_area",
    metric == "contig_mn" ~ "Contiguity_index"
  ))



saveRDS(urbanhmetrics1, "geo_data/lc_res5.rds")
saveRDS(area_res5, "geo_data/lc_area_res5.rds")
saveRDS(ce_res5, "geo_data/lc_condent_res5.rds")


#res3
hex_area_res3<-lapply(seq_along(res3_files), function(x){
  raster_read<-paste(directory_lc_res3, res3_files[[x]], sep = "/")
  raster_calc<-terra::rast(raster_read)
  calcs<-tryCatch(
    {calculate_lsm(raster_calc, what = c("lsm_l_ta"))},
    error = function(err) {return("error")}
  )
  
  rm(raster_calc)
  
  if(length(calcs)==1){ return(calcs)}else{
    calcs$res =str_extract(res3_files[[x]], "[^_]+")
    return(calcs)
  }
})

conditonal_entropy_res3<-lapply(seq_along(res3_files), function(x){
  raster_read<-paste(directory_lc_res3, res3_files[[x]], sep = "/")
  raster_calc<-terra::rast(raster_read)
  calcs<-tryCatch(
    {calculate_lsm(raster_calc, what = c("lsm_l_condent"))},
    error = function(err) {return("error")}
  )
  
  rm(raster_calc)
  
  if(length(calcs)==1){ return(calcs)}else{
    calcs$res =str_extract(res3_files[[x]], "[^_]+")
    return(calcs)
  }
})

hex_metrics_res3<-lapply(seq_along(res3_files), function(x){
  raster_read<-paste(directory_lc_res3, res3_files[[x]], sep = "/")
  raster_calc<-terra::rast(raster_read)
  plot(raster_calc)
  calcs<-tryCatch(
    {calculate_lsm(raster_calc, what = c("lsm_c_clumpy","lsm_c_ca",
                                         "lsm_c_contig_mn"))},
    
    error = function(err) {return("error")}
  )
  
  rm(raster_calc)
  
  if(length(calcs)==1){ return(calcs)}else{
    
    calcs$res = str_extract(res3_files[[x]], "[^_]+")
    return(calcs)
  }
  
})
remove<-which(hex_metrics_res3 == "error")
if(length(remove) == 0){
  conditonal_entropy_res3<-conditonal_entropy_res3
  hex_area_res3<-hex_area_res3
  hex_metrics_res3<-hex_metrics_res3
}else{
  conditonal_entropy_res3<-conditonal_entropy_res3[-remove]
  hex_area_res3<-hex_area_res3[-remove]
  hex_metrics_res3<-hex_metrics_res3[-remove]
}



ce_res3<-do.call(rbind, conditonal_entropy_res3) %>%
  dplyr::select(res, value) %>%
  rename(condent = value)
area_res3<-do.call(rbind, hex_area_res3) %>%
  dplyr::select(res, value)
metrics_res3<-do.call(rbind, hex_metrics_res3) %>%
  dplyr::select(res,class,metric, value)


urbanhmetrics1<-left_join(metrics_res3, ce_res3, by = "res") %>%
  mutate(metric = case_when(
    metric == "clumpy" ~ "Clumpiness_index",
    metric == "ca" ~ "Total_class_area",
    metric == "contig_mn" ~ "Contiguity_index"
  ))



saveRDS(urbanhmetrics1, "geo_data/lc_res3.rds")
saveRDS(area_res3, "geo_data/lc_area_res3.rds")
saveRDS(ce_res3, "geo_data/lc_condent_res3.rds")

