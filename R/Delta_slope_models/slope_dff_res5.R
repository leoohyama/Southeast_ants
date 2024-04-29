library(tidyverse)
library(sf)
library(h3jsr)
library(ape)


#sensitivity assessment for sample size in
#pgls

full_senseres4<-readRDS("pgls_results/FBD_crown_resampled_native_concensus_res5_sensitivity_data.rds")
full_senseres4_nn<-readRDS("pgls_results/FBD_crown_resampled_nonnative_concensus_res5_sensitivity_data.rds")

full_senseres4 <- drop_na(full_senseres4)
full_senseres4_nn <- drop_na(full_senseres4_nn)


ggplot(data = full_senseres4) +
  geom_point(aes(x = model_used_no_sig, y = mean_slope_QWD)) +
  facet_wrap(~sample_size)


ggplot(data = full_senseres4) +
  geom_point(aes(x = model_used_sig, y = mean_slope_QWD_sig)) +
  facet_wrap(~sample_size)


ggplot(data = full_senseres4_nn) +
  geom_point(aes(x = model_used_no_sig, y = mean_slope_QWD)) +
  facet_wrap(~sample_size)


ggplot(data = full_senseres4_nn) +
  geom_point(aes(x = model_used_sig, y = mean_slope_QWD_sig)) +
  facet_wrap(~sample_size)


#assess sample size effects

#read in landscapemetrics data
metrics_res4<-readRDS("geo_data/lc_res5.rds")
metrics_area_res4<-readRDS("geo_data/lc_area_res5.rds")
metrics_condent_res4<-readRDS("geo_data/lc_condent_res5.rds")


metrics_res4<-metrics_res4 %>%
  dplyr::select(class, metric, value, res) %>%
  filter(class %in% c(1,2,3,4)) %>%
  mutate(class = case_when(
    class == 1 ~ "Developed",
    class == 2 ~ "Cropland",
    class == 3 ~ "Grass_Shrub",
    class == 4 ~ "Tree_Cover"
  )) %>%
  pivot_wider(id_cols = res, names_from = c("class", "metric"), values_from = value)



metrics_condent_res4 <- metrics_condent_res4 %>%
  dplyr::select(res, condent) %>%
  rename(condent = condent)

metrics_area_res4 <- metrics_area_res4 %>%
  dplyr::select(res, value) %>%
  rename(total_area = value) %>%
  mutate(total_area = total_area/100)

native_only<-purrr::reduce(list(full_senseres4,
                                metrics_area_res4,
                                metrics_condent_res4, 
                                metrics_res4), 
                           dplyr::left_join, by = 'res')

#you need to separate by sample size
native_only<-native_only %>%
  group_by(sample_size) %>%
  nest()

native_nonnative<-purrr::reduce(list(full_senseres4_nn,
                                     metrics_area_res4,
                                     metrics_condent_res4, 
                                     metrics_res4), 
                                dplyr::left_join, by = 'res')

both_nn<-native_nonnative %>%
  group_by(sample_size) %>%
  nest()


intersecting_hex_20<-intersect(native_only$data[[1]]$res, both_nn$data[[1]]$res)
intersecting_hex_25<-intersect(native_only$data[[2]]$res, both_nn$data[[2]]$res)
intersecting_hex_30<-intersect(native_only$data[[3]]$res, both_nn$data[[3]]$res)

#now filter out by that list
list_intersect<-list(intersecting_hex_20,intersecting_hex_25,intersecting_hex_30)



native_cc<-map2(list_intersect, native_only$data, function(x,y){
  base<-y %>%
    select(c(res,model_used_no_sig,mean_r2,mean_intercept,mean_slope_QWD)) %>%
    filter(res %in% x) 
  centroid<-as_tibble(st_coordinates(st_centroid(h3jsr::cell_to_polygon(base$res))))
  base<-cbind(base, centroid)
  return(base)
})

native_nonnative_cc<-map2(list_intersect, both_nn$data, function(x,y){
  base<-y %>%
    select(c(res,model_used_no_sig,mean_r2,mean_intercept,mean_slope_QWD, exotic_perc_avg)) %>%
    filter(res %in% x) 
  centroid<-as_tibble(st_coordinates(st_centroid(h3jsr::cell_to_polygon(base$res))))
  base<-cbind(base, centroid)
  return(base)
})

#calculate basic difference
basic_diff<-map2(native_cc, native_nonnative_cc, function(x,y){
  x  = x %>% mutate(type = "native_only")
  y  = y %>% mutate(type = "both") %>%
    dplyr::select(-exotic_perc_avg)
  return(rbind(x,y))
})

ggplot(data = basic_diff[[1]]) +
  geom_boxplot(aes(x = type, y =mean_slope_QWD ))

ggplot(data = basic_diff[[2]]) +
  geom_boxplot(aes(x = type, y =mean_slope_QWD ))

ggplot(data = basic_diff[[3]]) +
  geom_boxplot(aes(x = type, y =mean_slope_QWD ))



differences<-map2(native_cc, native_nonnative_cc, function(x,y){
  diff_data<-left_join(x, y, by ="res")
  slope_diff = diff_data$mean_slope_QWD.y - diff_data$mean_slope_QWD.x
  r2_diff = diff_data$mean_r2.y - diff_data$mean_r2.x
  int_diff = diff_data$mean_intercept.y - diff_data$mean_intercept.x
  
  diff<-data.frame(res = diff_data$res, 
                   model_used_native = diff_data$model_used_no_sig.x,
                   model_used_both = diff_data$model_used_no_sig.y,
                   exotic_perc_avg = diff_data$exotic_perc_avg,
                   slope_diff = slope_diff,
                   r2_diff = r2_diff,
                   int_diff = int_diff,
                   X =diff_data$X.y, Y = diff_data$Y.y)
  diff<-purrr::reduce(list(diff,
                           metrics_area_res4,
                           metrics_condent_res4, 
                           metrics_res4), 
                      dplyr::left_join, by = 'res')
  return(diff)
  
})

saveRDS(differences, "difference_data/differences_res5.rds")


ggplot(differences[[1]]) +
  geom_point(aes(x = slope_diff, y= exotic_perc_avg*100)) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope difference", y = "Nonnative ant percent") +
  theme_light()



#sdmtmb part
library(sdmTMB)
library(exactextractr)
source("r/r2_sdmtmb_function.R")

#load climate
atr<-terra::rast("climate/wc2.0_bio_10m_07.tif")
mat<-terra::rast("climate/wc2.0_bio_10m_01.tif")

sample_20<-differences[[1]]
sample_25<-differences[[2]]
sample_30<-differences[[3]]

sample_20<-sample_20 %>%
  drop_na(Developed_Contiguity_index, Tree_Cover_Contiguity_index, 
          Tree_Cover_Clumpiness_index,Developed_Clumpiness_index) 
sample_25<-sample_25 %>%
  drop_na(Developed_Contiguity_index, Tree_Cover_Contiguity_index, 
          Tree_Cover_Clumpiness_index,Developed_Clumpiness_index) 
sample_30<-sample_30 %>%
  drop_na(Developed_Contiguity_index, Tree_Cover_Contiguity_index, 
          Tree_Cover_Clumpiness_index,Developed_Clumpiness_index) 


data_model<-sample_20
data_model25<-sample_25
data_model30<-sample_30

model_data<-list(data_model,data_model25,data_model30)

#add climate data to native data

model_data<-lapply(model_data, function(x){
  diff_poly<-cell_to_polygon(x$res) %>%
    st_as_sf()
  
  MAT <- cbind(diff_poly, exact_extract(mat, diff_poly, c('mean')))
  colnames(MAT)[1]<-"MAT"
  MAT = MAT %>% st_drop_geometry()
  ATR <- cbind(diff_poly, exact_extract(atr, diff_poly, c('mean')))
  colnames(ATR)[1]<-"ATR"
  ATR = ATR %>% st_drop_geometry()
  x<-cbind(x, data.frame(MAT =MAT, ATR= ATR))
  x$avg_model_used = (x$model_used_both+ x$model_used_native)/2
  return(x)
  
  
})

saveRDS(model_data, "difference_data/differences_res5.rds")


#check spatial auto for slope diff
lapply(model_data, function(x){
  inv_dists <- as.matrix(dist(x[,c("X","Y")]))
  diag(inv_dists) <- 0
  Moran.I(x$slope_diff, inv_dists)
})

hist(model_data[[1]]$slope_diff)
hist(model_data[[2]]$slope_diff)
hist(model_data[[3]]$slope_diff)


data_model<-model_data[[1]]
##now assess exotic percentage against slope sample size 20
meshe <- make_mesh(data_model, xy_cols = c("X", "Y"), n_knots = 60)

exotic_model <- sdmTMB(
  slope_diff ~  scale(avg_model_used) +
    scale(exotic_perc_avg),
  data = data_model,
  mesh = meshe,
  family = gaussian(link = "identity"),
  spatial = "on",
  silent = FALSE # see progress
)
sanity(exotic_model)

summary(exotic_model)
tidy(exotic_model,conf.int = TRUE,
     conf.level = 0.95)
summary(exotic_model$sd_report, select = "fixed", p.value = TRUE)
r2.sdmTMB(exotic_model)

#check residuals
s_gamma <- simulate(exotic_model, nsim = 500)
pred_fixed <- predict(exotic_model)
m1_resid <- DHARMa::createDHARMa(
  simulatedResponse = s_gamma,
  observedResponse = data_model$slope_diff,
  fittedPredictedResponse = pred_fixed$est
)
plot(m1_resid)


# significant effects of space found
Moran.I(data_model$r2_diff, inv_dists)
hist(data_model$r2_diff)

exotic_r2<-glm(r2_diff~scale(avg_model_used) +
                 scale(exotic_perc_avg), data = data_model)
summary(exotic_r2)

#check spatial auto for int
lapply(model_data, function(x){
  inv_dists <- as.matrix(dist(x[,c("X","Y")]))
  diag(inv_dists) <- 0
  Moran.I(x$int_diff, inv_dists)
})


exotic_model_int <- sdmTMB(
  int_diff ~  scale(avg_model_used) +
    scale(exotic_perc_avg),
  data = data_model,
  mesh = meshe,
  family = gaussian(link = "identity"),
  spatial = "on",
  silent = FALSE # see progress
)
sanity(exotic_model_int)
tidy(exotic_model_int,conf.int = TRUE,
     conf.level = 0.95)
summary(exotic_model_int$sd_report, select = "fixed", p.value = TRUE)
r2.sdmTMB(exotic_model_int)

#check spatial auto for r2
lapply(model_data, function(x){
  inv_dists <- as.matrix(dist(x[,c("X","Y")]))
  diag(inv_dists) <- 0
  Moran.I(x$r2_diff, inv_dists)
})
#no spatial autocorrelation in r2 differences

hist(data_model$r2_diff)

exotic_model_r2 <- glm(
  r2_diff ~  scale(avg_model_used) +
    scale(exotic_perc_avg),
  data = data_model,
  family = gaussian(link = "identity")
)
summary(exotic_model_r2)
confint(exotic_model_r2)
MuMIn::r.squaredGLMM(exotic_model_r2)



#now sample size 25

data_model<-model_data[[2]]
##now assess exotic percentage against slope sample size 20
meshe <- make_mesh(data_model, xy_cols = c("X", "Y"), n_knots = 45,type = "kmeans")

exotic_model <- sdmTMB(
  slope_diff ~  scale(avg_model_used) +
    scale(exotic_perc_avg),
  data = data_model,
  mesh = meshe,
  family = gaussian(link = "identity"),
  spatial = "on",
  silent = FALSE # see progress
)
sanity(exotic_model)
summary(exotic_model)
tidy(exotic_model)
tidy(exotic_model,conf.int = TRUE,
     conf.level = 0.95)
summary(exotic_model$sd_report, select = "fixed", p.value = TRUE)
r2.sdmTMB(exotic_model)

#check residuals
s_gamma <- simulate(exotic_model, nsim = 500)
pred_fixed <- predict(exotic_model)
m1_resid <- DHARMa::createDHARMa(
  simulatedResponse = s_gamma,
  observedResponse = data_model$slope_diff,
  fittedPredictedResponse = pred_fixed$est
)
plot(m1_resid)

#check spatial auto for int
lapply(model_data, function(x){
  inv_dists <- as.matrix(dist(x[,c("X","Y")]))
  diag(inv_dists) <- 0
  Moran.I(x$int_diff, inv_dists)
})


exotic_model_int <- sdmTMB(
  int_diff ~  scale(avg_model_used) +
    scale(exotic_perc_avg),
  data = data_model,
  mesh = meshe,
  family = gaussian(link = "identity"),
  spatial = "on",
  silent = FALSE # see progress
)
sanity(exotic_model_int)
tidy(exotic_model_int,conf.int = TRUE,
     conf.level = 0.95)
summary(exotic_model_int$sd_report, select = "fixed", p.value = TRUE)
r2.sdmTMB(exotic_model_int)

#check spatial auto for r2
lapply(model_data, function(x){
  inv_dists <- as.matrix(dist(x[,c("X","Y")]))
  diag(inv_dists) <- 0
  Moran.I(x$r2_diff, inv_dists)
})
#no spatial autocorrelation in r2 differences

hist(data_model$r2_diff)

exotic_model_r2 <- glm(
  r2_diff ~  scale(avg_model_used) +
    scale(exotic_perc_avg),
  data = data_model,
  family = gaussian(link = "identity")
)
summary(exotic_model_r2)
confint(exotic_model_r2)


MuMIn::r.squaredGLMM(exotic_model_r2)




#now sample size 30

data_model<-model_data[[3]]
##now assess exotic percentage against slope sample size 20
meshe <- make_mesh(data_model, xy_cols = c("X", "Y"), n_knots = 40,type = "kmeans")

exotic_model <- sdmTMB(
  slope_diff ~  scale(avg_model_used) +
    scale(exotic_perc_avg),
  data = data_model,
  mesh = meshe,
  family = gaussian(link = "identity"),
  spatial = "on",
  silent = FALSE # see progress
)
sanity(exotic_model)

summary(exotic_model)
tidy(exotic_model)
tidy(exotic_model,conf.int = TRUE,
     conf.level = 0.95)
summary(exotic_model$sd_report, select = "fixed", p.value = TRUE)
r2.sdmTMB(exotic_model)

#check residuals
s_gamma <- simulate(exotic_model, nsim = 500)
pred_fixed <- predict(exotic_model)
m1_resid <- DHARMa::createDHARMa(
  simulatedResponse = s_gamma,
  observedResponse = data_model$slope_diff,
  fittedPredictedResponse = pred_fixed$est
)
plot(m1_resid)

#check spatial auto for int
lapply(model_data, function(x){
  inv_dists <- as.matrix(dist(x[,c("X","Y")]))
  diag(inv_dists) <- 0
  Moran.I(x$int_diff, inv_dists)
})


exotic_model_int <- sdmTMB(
  int_diff ~  scale(avg_model_used) +
    scale(exotic_perc_avg),
  data = data_model,
  mesh = meshe,
  family = gaussian(link = "identity"),
  spatial = "on",
  silent = FALSE # see progress
)
sanity(exotic_model_int)
tidy(exotic_model_int,conf.int = TRUE,
     conf.level = 0.95)
summary(exotic_model_int$sd_report, select = "fixed", p.value = TRUE)
r2.sdmTMB(exotic_model_int)

#check residuals
s_gamma <- simulate(exotic_model_int, nsim = 500)
pred_fixed <- predict(exotic_model_int)
m1_resid <- DHARMa::createDHARMa(
  simulatedResponse = s_gamma,
  observedResponse = data_model$slope_diff,
  fittedPredictedResponse = pred_fixed$est
)
plot(m1_resid)

#check spatial auto for r2
lapply(model_data, function(x){
  inv_dists <- as.matrix(dist(x[,c("X","Y")]))
  diag(inv_dists) <- 0
  Moran.I(x$r2_diff, inv_dists)
})
#no spatial autocorrelation in r2 differences

hist(data_model$r2_diff)

exotic_model_r2 <- glm(
  r2_diff ~  scale(avg_model_used) +
    scale(exotic_perc_avg),
  data = data_model,
  family = gaussian(link = "identity")
)
summary(exotic_model_r2)
confint(exotic_model_r2)
MuMIn::r.squaredGLMM(exotic_model_r2)


