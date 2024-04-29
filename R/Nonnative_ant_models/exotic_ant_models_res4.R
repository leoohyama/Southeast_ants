library(tidyverse)
library(sf)
library(h3jsr)
library(exactextractr)
library(ape)
#sensitivity assessment for sample size in
#pgls

full_senseres4<-readRDS("pgls_results/FBD_crown_resampled_native_concensus_res4_sensitivity_data.rds")
full_senseres4_nn<-readRDS("pgls_results/FBD_crown_resampled_nonnative_concensus_res4_sensitivity_data.rds")

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
metrics_res4<-readRDS("geo_data/lc_res4.rds")
metrics_area_res4<-readRDS("geo_data/lc_area_res4.rds")
metrics_condent_res4<-readRDS("geo_data/lc_condent_res4.rds")

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
  dplyr::select(res, condent) 

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
saveRDS(differences, "difference_data/differences_res4.rds")




#sdmtmb part
library(sdmTMB)
source("r/r2_sdmtmb_function.R")

#load climate
atr<-terra::rast("climate/wc2.0_bio_10m_07.tif")
mat<-terra::rast("climate/wc2.0_bio_10m_01.tif")

sample_20<-differences[[1]]
sample_25<-differences[[2]]
sample_30<-differences[[3]]

sample_20<-sample_20 %>%
  drop_na(Developed_Contiguity_index, Tree_Cover_Contiguity_index, 
          Tree_Cover_Clumpiness_index,Developed_Clumpiness_index) %>%
  filter(!exotic_perc_avg == 0)
sample_25<-sample_25 %>%
  drop_na(Developed_Contiguity_index, Tree_Cover_Contiguity_index, 
          Tree_Cover_Clumpiness_index,Developed_Clumpiness_index) %>%
  filter(!exotic_perc_avg == 0)
sample_30<-sample_30 %>%
  drop_na(Developed_Contiguity_index, Tree_Cover_Contiguity_index, 
          Tree_Cover_Clumpiness_index,Developed_Clumpiness_index) %>%
  filter(!exotic_perc_avg == 0)

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
#check spatial auto for slope diff
lapply(model_data, function(x){
  inv_dists <- as.matrix(dist(x[,c("X","Y")]))
  diag(inv_dists) <- 0
  Moran.I(x$exotic_perc_avg, inv_dists)
})

hist(model_data[[1]]$exotic_perc_avg)
hist(model_data[[2]]$exotic_perc_avg)
hist(model_data[[3]]$exotic_perc_avg)


data_model<-model_data[[1]]
mesh <- make_mesh(data_model, xy_cols = c("X", "Y"), n_knots = 60, type = "kmeans")
plot(mesh)

hist(data_model$exotic_perc_avg)


dev_model_cont_beta <- sdmTMB(
  exotic_perc_avg ~ scale(model_used_both) +
    scale(ATR) +
    scale(Developed_Contiguity_index),
  data = data_model,
  mesh = mesh,
  family = gaussian(link = "identity"),
  spatial = "on",
  silent = FALSE # see progress,
  
)


dev_model_cont_beta_int <- sdmTMB(
  (exotic_perc_avg) ~ scale(model_used_both) +
    scale(ATR) *
    scale(Developed_Contiguity_index),
  data = data_model,
  mesh = mesh,
  family = Beta(link = "logit"),
  control = sdmTMBcontrol(),
  spatial = "on",
  silent = FALSE # see progress
)

dev_model_clump_beta <- sdmTMB(
  exotic_perc_avg ~  scale(model_used_both) +
    scale(ATR) +
    scale(Developed_Clumpiness_index),
  data = data_model,
  mesh = mesh,
  family = Beta(link = "logit"),
  spatial = "on",
  silent = FALSE # see progress
)

dev_model_clump_beta_int<- sdmTMB(
  (exotic_perc_avg) ~  scale(model_used_both) +
    scale(ATR) *
    scale(Developed_Clumpiness_index),
  data = data_model,
  mesh = mesh,
  family = Beta(link = "logit"),
  spatial = "on",
  silent = FALSE # see progress
)


null_mod <- sdmTMB(
  exotic_perc_avg ~ 
    1,
  data = data_model,
  mesh = mesh,
  family = Beta(link = "logit"),
  spatial = "on",
  silent = FALSE # see progress
)
library(bbmle)
AICtab(dev_model_clump_beta_int,dev_model_clump_beta,dev_model_cont_beta,
       dev_model_cont_beta_int, null_mod, weights = T, base =T)



bestmodel<-dev_model_clump_beta
summary(bestmodel)
summary(bestmodel$sd_report, select = "fixed", p.value = TRUE)
tidy(bestmodel, conf.int = T, conf.level = 0.95)

library(DHARMa)
#check residuals
s_gamma <- simulate(bestmodel, nsim = 500)
pred_fixed <- bestmodel$family$linkinv(predict(bestmodel)$est_non_rf)
m1_resid <- DHARMa::createDHARMa(
  simulatedResponse = s_gamma,
  observedResponse = data_model$exotic_perc_avg,
  fittedPredictedResponse = pred_fixed
)

plot(m1_resid)
inv_dists <- as.matrix(dist(data_model[,c("X","Y")]))
diag(inv_dists) <- 0
Moran.I(resid(bestmodel), inv_dists)


data_model<-model_data[[2]]
mesh <- make_mesh(data_model, xy_cols = c("X", "Y"), n_knots = 50, type = "kmeans")
plot(mesh)


dev_model_cont_beta <- sdmTMB(
  exotic_perc_avg ~ scale(model_used_both) +
    scale(ATR) +
    scale(Developed_Contiguity_index),
  data = data_model,
  mesh = mesh,
  family = Beta(link = "logit"),
  spatial = "on",
  silent = FALSE # see progress,
  
)


dev_model_cont_beta_int <- sdmTMB(
  (exotic_perc_avg) ~ scale(model_used_both) +
    scale(ATR) *
    scale(Developed_Contiguity_index),
  data = data_model,
  mesh = mesh,
  family = Beta(link = "logit"),
  control = sdmTMBcontrol(),
  spatial = "on",
  silent = FALSE # see progress
)

dev_model_clump_beta <- sdmTMB(
  exotic_perc_avg ~  scale(model_used_both) +
    scale(ATR) +
    scale(Developed_Clumpiness_index),
  data = data_model,
  mesh = mesh,
  family = Beta(link = "logit"),
  spatial = "on",
  silent = FALSE # see progress
)

dev_model_clump_beta_int<- sdmTMB(
  (exotic_perc_avg) ~  scale(model_used_both) +
    scale(ATR) *
    scale(Developed_Clumpiness_index),
  data = data_model,
  mesh = mesh,
  family = Beta(link = "logit"),
  spatial = "on",
  silent = FALSE # see progress
)


null_mod <- sdmTMB(
  exotic_perc_avg ~ 
    1,
  data = data_model,
  mesh = mesh,
  family = gaussian(link = "identity"),
  spatial = "on",
  silent = FALSE # see progress
)
library(bbmle)
AICtab(dev_model_clump_beta_int,dev_model_clump_beta,dev_model_cont_beta,
       dev_model_cont_beta_int, null_mod, weights = T, base = T)


bestmodel<-dev_model_clump_beta

summary(bestmodel)
summary(bestmodel$sd_report, select = "fixed", p.value = TRUE)
tidy(bestmodel, conf.int = T, conf.level = 0.95)

library(DHARMa)
#check residuals
s_gamma <- simulate(bestmodel, nsim = 500)
pred_fixed <- bestmodel$family$linkinv(predict(bestmodel)$est_non_rf)
m1_resid <- DHARMa::createDHARMa(
  simulatedResponse = s_gamma,
  observedResponse = data_model$exotic_perc_avg,
  fittedPredictedResponse = pred_fixed
)

plot(m1_resid)
inv_dists <- as.matrix(dist(data_model[,c("X","Y")]))
diag(inv_dists) <- 0
Moran.I(resid(bestmodel), inv_dists)

data_model<-model_data[[3]]
mesh <- make_mesh(data_model, xy_cols = c("X", "Y"), n_knots = 40, type = "kmeans")
plot(mesh)


dev_model_cont_beta <- sdmTMB(
  exotic_perc_avg ~ scale(avg_model_used) +
    scale(ATR) +
    scale(Developed_Contiguity_index),
  data = data_model,
  mesh = mesh,
  family = Beta(link = "logit"),
  spatial = "on",
  silent = FALSE # see progress,
  
)


dev_model_cont_beta_int <- sdmTMB(
  (exotic_perc_avg) ~ scale(avg_model_used) +
    scale(ATR) *
    scale(Developed_Contiguity_index),
  data = data_model,
  mesh = mesh,
  family = Beta(link = "logit"),
  control = sdmTMBcontrol(),
  spatial = "on",
  silent = FALSE # see progress
)

dev_model_clump_beta <- sdmTMB(
  exotic_perc_avg ~  scale(avg_model_used) +
    scale(ATR) +
    scale(Developed_Clumpiness_index),
  data = data_model,
  mesh = mesh,
  family = Beta(link = "logit"),
  spatial = "on",
  silent = FALSE # see progress
)

dev_model_clump_beta_int<- sdmTMB(
  (exotic_perc_avg) ~  scale(avg_model_used) +
    scale(ATR) *
    scale(Developed_Clumpiness_index),
  data = data_model,
  mesh = mesh,
  family = Beta(link = "logit"),
  spatial = "on",
  silent = FALSE # see progress
)


null_mod <- sdmTMB(
  exotic_perc_avg ~ 
    1,
  data = data_model,
  mesh = mesh,
  family = gaussian(link = "identity"),
  spatial = "on",
  silent = FALSE # see progress
)
library(bbmle)
AICtab(dev_model_clump_beta_int,dev_model_clump_beta,dev_model_cont_beta,
       dev_model_cont_beta_int, null_mod, weights = T, base =T)


bestmodel<-dev_model_clump_beta
summary(bestmodel)
sanity(bestmodel)
summary(bestmodel$sd_report, select = "fixed", p.value = TRUE)
tidy(bestmodel, conf.int = T, conf.level = 0.95)

library(DHARMa)
#check residuals
s_gamma <- simulate(bestmodel, nsim = 500)
pred_fixed <- bestmodel$family$linkinv(predict(bestmodel)$est_non_rf)
m1_resid <- DHARMa::createDHARMa(
  simulatedResponse = s_gamma,
  observedResponse = data_model$exotic_perc_avg,
  fittedPredictedResponse = pred_fixed
)

plot(m1_resid)
inv_dists <- as.matrix(dist(data_model[,c("X","Y")]))
diag(inv_dists) <- 0
Moran.I(resid(bestmodel), inv_dists)


