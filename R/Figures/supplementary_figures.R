library(tidyverse)
library(patchwork)
#sensitivity assessment for sample size in
#pgls

full_senseres4<-readRDS("pgls_results/FBD_crown_resampled_native_concensus_res4_sensitivity_data.rds")
full_senseres4_nn<-readRDS("pgls_results/FBD_crown_resampled_nonnative_concensus_res4_sensitivity_data.rds")

full_senseres4 <- drop_na(full_senseres4) %>%
  mutate(resolution = "~1,770km2")
full_senseres4_nn <- drop_na(full_senseres4_nn)%>%
  mutate(resolution = "~1,770km2")

full_senseres5<-readRDS("pgls_results/FBD_crown_resampled_native_concensus_res5_sensitivity_data.rds")
full_senseres5_nn<-readRDS("pgls_results/FBD_crown_resampled_nonnative_concensus_res5_sensitivity_data.rds")

full_senseres5 <- drop_na(full_senseres5)%>%
  mutate(resolution = "~252km2")
full_senseres5_nn <- drop_na(full_senseres5_nn)%>%
  mutate(resolution = "~252km2")


full_senseres3<-readRDS("pgls_results/FBD_crown_resampled_native_concensus_res3_sensitivity_data.rds")
full_senseres3_nn<-readRDS("pgls_results/FBD_crown_resampled_nonnative_concensus_res3_sensitivity_data.rds")

full_senseres3 <- drop_na(full_senseres3)%>%
  mutate(resolution = "~12,330km2")
full_senseres3_nn <- drop_na(full_senseres3_nn)%>%
  mutate(resolution = "~12,330km2")



native_only<-rbind(full_senseres4,full_senseres5,full_senseres3)
nativeandnon<-rbind(full_senseres4_nn,full_senseres5_nn,full_senseres3_nn)





native_only_QWD<-ggplot(data = native_only) +
  geom_point(aes(x = model_used_no_sig, y = mean_slope_QWD)) +
  facet_grid(vars(sample_size), vars(resolution)) +
  labs(x = "No. of PGLS models to calculate parameter",
       y = "Average QWD slope coefficient",title = "Native Species only data") +
  theme_bw()


nonnative_native_QWD<-ggplot(data = nativeandnon) +
  geom_point(aes(x = model_used_no_sig, y = mean_slope_QWD)) +
  facet_grid(vars(sample_size), vars(resolution)) +
  labs(x = "No. of PGLS models to calculate parameter",
       y = "Average QWD slope coefficient",title = "Native Species and Nonnative Species") +
  theme_bw()

native_only_QWD +nonnative_native_QWD


native_only_r2<-ggplot(data = native_only) +
  geom_point(aes(x = model_used_no_sig, y = mean_r2)) +
  facet_grid(vars(sample_size), vars(resolution)) +
  labs(x = "No. of PGLS models to calculate parameter",
       y = "Average R2",title = "Native Species only data") +
  theme_bw()


nonnative_native_r2<-ggplot(data = nativeandnon) +
  geom_point(aes(x = model_used_no_sig, y = mean_r2)) +
  facet_grid(vars(sample_size), vars(resolution)) +
  labs(x = "No. of PGLS models to calculate parameter",
       y = "Average R2 ",title = "Native Species and Nonnative Species") +
  theme_bw()

native_only_r2 +nonnative_native_r2


native_only_int<-ggplot(data = native_only) +
  geom_point(aes(x = model_used_no_sig, y = mean_intercept)) +
  facet_grid(vars(sample_size), vars(resolution)) +
  labs(x = "No. of PGLS models to calculate parameter",
       y = "Average intercept",title = "Native Species only data") +
  theme_bw()


nonnative_native_int<-ggplot(data = nativeandnon) +
  geom_point(aes(x = model_used_no_sig, y = mean_intercept)) +
  facet_grid(vars(sample_size), vars(resolution)) +
  labs(x = "No. of PGLS models to calculate parameter",
       y = "Average intercept",title = "Native Species and Nonnative Species") +
  theme_bw()

native_only_int +nonnative_native_int

native_only1<-native_only %>% mutate(type = "Native Species only")
nativeandnon1<-nativeandnon %>% mutate(type = "Native Species + NonnativeSpecies") %>%
  dplyr::select(-c(exotic_perc_avg_sig,exotic_perc_avg))

total_data<-rbind(native_only1, nativeandnon1)


slopes_mean<-ggplot(data = total_data) +
  geom_boxplot(aes(x = as.factor(sample_size), y = mean_slope_QWD)) +
  facet_grid(vars(resolution), vars(type)) +
  labs(x = "Sample Size PGLS",
       y = "Average QWD slope coefficient") +
  theme_bw()


r2_mean<-ggplot(data = total_data) +
  geom_boxplot(aes(x = as.factor(sample_size), y = mean_r2)) +
  facet_grid(vars(resolution), vars(type)) +
  labs(x = "Sample Size PGLS",
       y = "Average R2") +
  theme_bw()


int_mean<-ggplot(data = total_data) +
  geom_boxplot(aes(x = as.factor(sample_size), y = mean_intercept)) +
  facet_grid(vars(resolution), vars(type)) +
  labs(x = "Sample Size PGLS",
       y = "Average Intercept") +
  theme_bw()



null_mod<-ggplot(data = total_data) +
  geom_boxplot(aes(x = as.factor(sample_size), y = n_null)) +
  facet_grid(vars(resolution), vars(type)) +
  labs(x = "Sample Size PGLS",
       y = "No. iterations (out of 100) when the null model is most plausible") +
  theme_bw()



#comparing PGLS parameters across different tree types
library(tidyverse)
#Native pgls @ res 3

FBD_crown_20_native_res3<-readRDS("pgls_results/FBD_crown_resampled_20_native_concensus_res3.rds")
FBD_crown_20_native_res3$Tree<-"FBD_crown"
FBD_stem_20_native_res3<-readRDS("pgls_results/FBD_stem_resampled_20_native_concensus_res3.rds")
FBD_stem_20_native_res3$Tree<-"FBD_stem"
NC_crown_20_native_res3<-readRDS("pgls_results/NC_crown_resampled_20_native_concensus_res3.rds")
NC_crown_20_native_res3$Tree<-"NC_crown"
NC_stem_20_native_res3<-readRDS("pgls_results/NC_stem_resampled_20_native_concensus_res3.rds")
NC_stem_20_native_res3$Tree<-"NC_stem"
res3_native_20<-rbind(FBD_crown_20_native_res3, FBD_stem_20_native_res3,
                      NC_crown_20_native_res3, NC_stem_20_native_res3)
res3_native_20$sample_size<-20


FBD_crown_25_native_res3<-readRDS("pgls_results/FBD_crown_resampled_25_native_concensus_res3.rds")
FBD_crown_25_native_res3$Tree<-"FBD_crown"
FBD_stem_25_native_res3<-readRDS("pgls_results/FBD_stem_resampled_25_native_concensus_res3.rds")
FBD_stem_25_native_res3$Tree<-"FBD_stem"
NC_crown_25_native_res3<-readRDS("pgls_results/NC_crown_resampled_25_native_concensus_res3.rds")
NC_crown_25_native_res3$Tree<-"NC_crown"
NC_stem_25_native_res3<-readRDS("pgls_results/NC_stem_resampled_25_native_concensus_res3.rds")
NC_stem_25_native_res3$Tree<-"NC_stem"

res3_native_25<-rbind(FBD_crown_25_native_res3, FBD_stem_25_native_res3,
                      NC_crown_25_native_res3, NC_stem_25_native_res3)
res3_native_25$sample_size<-25



FBD_crown_30_native_res3<-readRDS("pgls_results/FBD_crown_resampled_30_native_concensus_res3.rds")
FBD_crown_30_native_res3$Tree<-"FBD_crown"
FBD_stem_30_native_res3<-readRDS("pgls_results/FBD_stem_resampled_30_native_concensus_res3.rds")
FBD_stem_30_native_res3$Tree<-"FBD_stem"
NC_crown_30_native_res3<-readRDS("pgls_results/NC_crown_resampled_30_native_concensus_res3.rds")
NC_crown_30_native_res3$Tree<-"NC_crown"
NC_stem_30_native_res3<-readRDS("pgls_results/NC_stem_resampled_30_native_concensus_res3.rds")
NC_stem_30_native_res3$Tree<-"NC_stem"

res3_native_30<-rbind(FBD_crown_30_native_res3, FBD_stem_30_native_res3,
                      NC_crown_30_native_res3, NC_stem_30_native_res3)
res3_native_30$sample_size<-30

res3_all_samples_native<-rbind(res3_native_20,res3_native_25,res3_native_30)
res3_all_samples_native$status<-"native_species_only"


#Native + nonnative pgls @ res 3


FBD_crown_20_nonnative_res3<-readRDS("pgls_results/FBD_crown_resampled_20_nonnative_concensus_res3.rds")
FBD_crown_20_nonnative_res3<-FBD_crown_20_nonnative_res3[-14]
FBD_crown_20_nonnative_res3$Tree<-"FBD_crown"
FBD_stem_20_nonnative_res3<-readRDS("pgls_results/FBD_stem_resampled_20_nonnative_concensus_res3.rds")
FBD_stem_20_nonnative_res3<-FBD_stem_20_nonnative_res3[-14]
FBD_stem_20_nonnative_res3$Tree<-"FBD_stem"
NC_crown_20_nonnative_res3<-readRDS("pgls_results/NC_crown_resampled_20_nonnative_concensus_res3.rds")
NC_crown_20_nonnative_res3<-NC_crown_20_nonnative_res3[-14]
NC_crown_20_nonnative_res3$Tree<-"NC_crown"
NC_stem_20_nonnative_res3<-readRDS("pgls_results/NC_stem_resampled_20_nonnative_concensus_res3.rds")
NC_stem_20_nonnative_res3<-NC_stem_20_nonnative_res3[-14]
NC_stem_20_nonnative_res3$Tree<-"NC_stem"
res3_nonnative_20<-rbind(FBD_crown_20_nonnative_res3, FBD_stem_20_nonnative_res3,
                         NC_crown_20_nonnative_res3, NC_stem_20_nonnative_res3)
res3_nonnative_20$sample_size<-20


FBD_crown_25_nonnative_res3<-readRDS("pgls_results/FBD_crown_resampled_25_nonnative_concensus_res3.rds")
FBD_crown_25_nonnative_res3<-FBD_crown_25_nonnative_res3[-14]
FBD_crown_25_nonnative_res3$Tree<-"FBD_crown"
FBD_stem_25_nonnative_res3<-readRDS("pgls_results/FBD_stem_resampled_25_nonnative_concensus_res3.rds")
FBD_stem_25_nonnative_res3<-FBD_stem_25_nonnative_res3[-14]
FBD_stem_25_nonnative_res3$Tree<-"FBD_stem"
NC_crown_25_nonnative_res3<-readRDS("pgls_results/NC_crown_resampled_25_nonnative_concensus_res3.rds")
NC_crown_25_nonnative_res3<-NC_crown_25_nonnative_res3[-14]
NC_crown_25_nonnative_res3$Tree<-"NC_crown"
NC_stem_25_nonnative_res3<-readRDS("pgls_results/NC_stem_resampled_25_nonnative_concensus_res3.rds")
NC_stem_25_nonnative_res3<-NC_stem_25_nonnative_res3[-14]
NC_stem_25_nonnative_res3$Tree<-"NC_stem"

res3_nonnative_25<-rbind(FBD_crown_25_nonnative_res3, FBD_stem_25_nonnative_res3,
                         NC_crown_25_nonnative_res3, NC_stem_25_nonnative_res3)
res3_nonnative_25$sample_size<-25



FBD_crown_30_nonnative_res3<-readRDS("pgls_results/FBD_crown_resampled_30_nonnative_concensus_res3.rds")
FBD_crown_30_nonnative_res3<-FBD_crown_30_nonnative_res3[-14]
FBD_crown_30_nonnative_res3$Tree<-"FBD_crown"
FBD_stem_30_nonnative_res3<-readRDS("pgls_results/FBD_stem_resampled_30_nonnative_concensus_res3.rds")
FBD_stem_30_nonnative_res3<-FBD_stem_30_nonnative_res3[-14]
FBD_stem_30_nonnative_res3$Tree<-"FBD_stem"
NC_crown_30_nonnative_res3<-readRDS("pgls_results/NC_crown_resampled_30_nonnative_concensus_res3.rds")
NC_crown_30_nonnative_res3<-NC_crown_30_nonnative_res3[-14]
NC_crown_30_nonnative_res3$Tree<-"NC_crown"
NC_stem_30_nonnative_res3<-readRDS("pgls_results/NC_stem_resampled_30_nonnative_concensus_res3.rds")
NC_stem_30_nonnative_res3<-NC_stem_30_nonnative_res3[-14]
NC_stem_30_nonnative_res3$Tree<-"NC_stem"

res3_nonnative_30<-rbind(FBD_crown_30_nonnative_res3, FBD_stem_30_nonnative_res3,
                         NC_crown_30_nonnative_res3, NC_stem_30_nonnative_res3)
res3_nonnative_30$sample_size<-30

res3_all_samples_nonnative<-rbind(res3_nonnative_20,res3_nonnative_25,res3_nonnative_30)
res3_all_samples_nonnative$status<-"native_nonnative_species"

res3_all<-rbind(res3_all_samples_nonnative, res3_all_samples_native)

res3_all<-res3_all %>% filter(!r2 == 0)

res3_all %>%
  mutate(slope_QWD = as.numeric(slope_QWD),
         sample_size = as.factor(sample_size)) %>%
  ggplot(.) +
  geom_boxplot(aes(x = sample_size, y = slope_QWD, fill = Tree)) +
  facet_wrap(~status)+
  labs(title = "@ 12,930km2")



#comparing PGLS parameters across different tree types
library(tidyverse)
#Native pgls @ res 4

FBD_crown_20_native_res4<-readRDS("pgls_results/FBD_crown_resampled_20_native_concensus_res4.rds")
FBD_crown_20_native_res4$Tree<-"FBD_crown"
FBD_stem_20_native_res4<-readRDS("pgls_results/FBD_stem_resampled_20_native_concensus_res4.rds")
FBD_stem_20_native_res4$Tree<-"FBD_stem"
NC_crown_20_native_res4<-readRDS("pgls_results/NC_crown_resampled_20_native_concensus_res4.rds")
NC_crown_20_native_res4$Tree<-"NC_crown"
NC_stem_20_native_res4<-readRDS("pgls_results/NC_stem_resampled_20_native_concensus_res4.rds")
NC_stem_20_native_res4$Tree<-"NC_stem"
res4_native_20<-rbind(FBD_crown_20_native_res4, FBD_stem_20_native_res4,
                      NC_crown_20_native_res4, NC_stem_20_native_res4)
res4_native_20$sample_size<-20


FBD_crown_25_native_res4<-readRDS("pgls_results/FBD_crown_resampled_25_native_concensus_res4.rds")
FBD_crown_25_native_res4$Tree<-"FBD_crown"
FBD_stem_25_native_res4<-readRDS("pgls_results/FBD_stem_resampled_25_native_concensus_res4.rds")
FBD_stem_25_native_res4$Tree<-"FBD_stem"
NC_crown_25_native_res4<-readRDS("pgls_results/NC_crown_resampled_25_native_concensus_res4.rds")
NC_crown_25_native_res4$Tree<-"NC_crown"
NC_stem_25_native_res4<-readRDS("pgls_results/NC_stem_resampled_25_native_concensus_res4.rds")
NC_stem_25_native_res4$Tree<-"NC_stem"

res4_native_25<-rbind(FBD_crown_25_native_res4, FBD_stem_25_native_res4,
                      NC_crown_25_native_res4, NC_stem_25_native_res4)
res4_native_25$sample_size<-25



FBD_crown_30_native_res4<-readRDS("pgls_results/FBD_crown_resampled_30_native_concensus_res4.rds")
FBD_crown_30_native_res4$Tree<-"FBD_crown"
FBD_stem_30_native_res4<-readRDS("pgls_results/FBD_stem_resampled_30_native_concensus_res4.rds")
FBD_stem_30_native_res4$Tree<-"FBD_stem"
NC_crown_30_native_res4<-readRDS("pgls_results/NC_crown_resampled_30_native_concensus_res4.rds")
NC_crown_30_native_res4$Tree<-"NC_crown"
NC_stem_30_native_res4<-readRDS("pgls_results/NC_stem_resampled_30_native_concensus_res4.rds")
NC_stem_30_native_res4$Tree<-"NC_stem"

res4_native_30<-rbind(FBD_crown_30_native_res4, FBD_stem_30_native_res4,
                      NC_crown_30_native_res4, NC_stem_30_native_res4)
res4_native_30$sample_size<-30

res4_all_samples_native<-rbind(res4_native_20,res4_native_25,res4_native_30)
res4_all_samples_native$status<-"native_species_only"


#Native + nonnative pgls @ res 4


FBD_crown_20_nonnative_res4<-readRDS("pgls_results/FBD_crown_resampled_20_nonnative_concensus_res4.rds")
FBD_crown_20_nonnative_res4<-FBD_crown_20_nonnative_res4[-14]
FBD_crown_20_nonnative_res4$Tree<-"FBD_crown"
FBD_stem_20_nonnative_res4<-readRDS("pgls_results/FBD_stem_resampled_20_nonnative_concensus_res4.rds")
FBD_stem_20_nonnative_res4<-FBD_stem_20_nonnative_res4[-14]
FBD_stem_20_nonnative_res4$Tree<-"FBD_stem"
NC_crown_20_nonnative_res4<-readRDS("pgls_results/NC_crown_resampled_20_nonnative_concensus_res4.rds")
NC_crown_20_nonnative_res4<-NC_crown_20_nonnative_res4[-14]
NC_crown_20_nonnative_res4$Tree<-"NC_crown"
NC_stem_20_nonnative_res4<-readRDS("pgls_results/NC_stem_resampled_20_nonnative_concensus_res4.rds")
NC_stem_20_nonnative_res4<-NC_stem_20_nonnative_res4[-14]
NC_stem_20_nonnative_res4$Tree<-"NC_stem"
res4_nonnative_20<-rbind(FBD_crown_20_nonnative_res4, FBD_stem_20_nonnative_res4,
                         NC_crown_20_nonnative_res4, NC_stem_20_nonnative_res4)
res4_nonnative_20$sample_size<-20


FBD_crown_25_nonnative_res4<-readRDS("pgls_results/FBD_crown_resampled_25_nonnative_concensus_res4.rds")
FBD_crown_25_nonnative_res4<-FBD_crown_25_nonnative_res4[-14]
FBD_crown_25_nonnative_res4$Tree<-"FBD_crown"
FBD_stem_25_nonnative_res4<-readRDS("pgls_results/FBD_stem_resampled_25_nonnative_concensus_res4.rds")
FBD_stem_25_nonnative_res4<-FBD_stem_25_nonnative_res4[-14]
FBD_stem_25_nonnative_res4$Tree<-"FBD_stem"
NC_crown_25_nonnative_res4<-readRDS("pgls_results/NC_crown_resampled_25_nonnative_concensus_res4.rds")
NC_crown_25_nonnative_res4<-NC_crown_25_nonnative_res4[-14]
NC_crown_25_nonnative_res4$Tree<-"NC_crown"
NC_stem_25_nonnative_res4<-readRDS("pgls_results/NC_stem_resampled_25_nonnative_concensus_res4.rds")
NC_stem_25_nonnative_res4<-NC_stem_25_nonnative_res4[-14]
NC_stem_25_nonnative_res4$Tree<-"NC_stem"

res4_nonnative_25<-rbind(FBD_crown_25_nonnative_res4, FBD_stem_25_nonnative_res4,
                         NC_crown_25_nonnative_res4, NC_stem_25_nonnative_res4)
res4_nonnative_25$sample_size<-25



FBD_crown_30_nonnative_res4<-readRDS("pgls_results/FBD_crown_resampled_30_nonnative_concensus_res4.rds")
FBD_crown_30_nonnative_res4<-FBD_crown_30_nonnative_res4[-14]
FBD_crown_30_nonnative_res4$Tree<-"FBD_crown"
FBD_stem_30_nonnative_res4<-readRDS("pgls_results/FBD_stem_resampled_30_nonnative_concensus_res4.rds")
FBD_stem_30_nonnative_res4<-FBD_stem_30_nonnative_res4[-14]
FBD_stem_30_nonnative_res4$Tree<-"FBD_stem"
NC_crown_30_nonnative_res4<-readRDS("pgls_results/NC_crown_resampled_30_nonnative_concensus_res4.rds")
NC_crown_30_nonnative_res4<-NC_crown_30_nonnative_res4[-14]
NC_crown_30_nonnative_res4$Tree<-"NC_crown"
NC_stem_30_nonnative_res4<-readRDS("pgls_results/NC_stem_resampled_30_nonnative_concensus_res4.rds")
NC_stem_30_nonnative_res4<-NC_stem_30_nonnative_res4[-14]
NC_stem_30_nonnative_res4$Tree<-"NC_stem"

res4_nonnative_30<-rbind(FBD_crown_30_nonnative_res4, FBD_stem_30_nonnative_res4,
                         NC_crown_30_nonnative_res4, NC_stem_30_nonnative_res4)
res4_nonnative_30$sample_size<-30

res4_all_samples_nonnative<-rbind(res4_nonnative_20,res4_nonnative_25,res4_nonnative_30)
res4_all_samples_nonnative$status<-"native_nonnative_species"

res4_all<-rbind(res4_all_samples_nonnative, res4_all_samples_native)

res4_all<-res4_all %>% filter(!r2 == 0)

res4_all %>%
  mutate(slope_QWD = as.numeric(slope_QWD),
         sample_size = as.factor(sample_size)) %>%
  ggplot(.) +
  geom_boxplot(aes(x = sample_size, y = slope_QWD, fill = Tree)) +
  facet_wrap(~status) +
  labs(title = "@ 1,770km2")


#comparing PGLS parameters across different tree types
library(tidyverse)
#Native pgls @ res 5

FBD_crown_20_native_res5<-readRDS("pgls_results/FBD_crown_resampled_20_native_concensus_res5.rds")
FBD_crown_20_native_res5$Tree<-"FBD_crown"
FBD_stem_20_native_res5<-readRDS("pgls_results/FBD_stem_resampled_20_native_concensus_res5.rds")
FBD_stem_20_native_res5$Tree<-"FBD_stem"
NC_crown_20_native_res5<-readRDS("pgls_results/NC_crown_resampled_20_native_concensus_res5.rds")
NC_crown_20_native_res5$Tree<-"NC_crown"
NC_stem_20_native_res5<-readRDS("pgls_results/NC_stem_resampled_20_native_concensus_res5.rds")
NC_stem_20_native_res5$Tree<-"NC_stem"
res5_native_20<-rbind(FBD_crown_20_native_res5, FBD_stem_20_native_res5,
                      NC_crown_20_native_res5, NC_stem_20_native_res5)
res5_native_20$sample_size<-20


FBD_crown_25_native_res5<-readRDS("pgls_results/FBD_crown_resampled_25_native_concensus_res5.rds")
FBD_crown_25_native_res5$Tree<-"FBD_crown"
FBD_stem_25_native_res5<-readRDS("pgls_results/FBD_stem_resampled_25_native_concensus_res5.rds")
FBD_stem_25_native_res5$Tree<-"FBD_stem"
NC_crown_25_native_res5<-readRDS("pgls_results/NC_crown_resampled_25_native_concensus_res5.rds")
NC_crown_25_native_res5$Tree<-"NC_crown"
NC_stem_25_native_res5<-readRDS("pgls_results/NC_stem_resampled_25_native_concensus_res5.rds")
NC_stem_25_native_res5$Tree<-"NC_stem"

res5_native_25<-rbind(FBD_crown_25_native_res5, FBD_stem_25_native_res5,
                      NC_crown_25_native_res5, NC_stem_25_native_res5)
res5_native_25$sample_size<-25



FBD_crown_30_native_res5<-readRDS("pgls_results/FBD_crown_resampled_30_native_concensus_res5.rds")
FBD_crown_30_native_res5$Tree<-"FBD_crown"
FBD_stem_30_native_res5<-readRDS("pgls_results/FBD_stem_resampled_30_native_concensus_res5.rds")
FBD_stem_30_native_res5$Tree<-"FBD_stem"
NC_crown_30_native_res5<-readRDS("pgls_results/NC_crown_resampled_30_native_concensus_res5.rds")
NC_crown_30_native_res5$Tree<-"NC_crown"
NC_stem_30_native_res5<-readRDS("pgls_results/NC_stem_resampled_30_native_concensus_res5.rds")
NC_stem_30_native_res5$Tree<-"NC_stem"

res5_native_30<-rbind(FBD_crown_30_native_res5, FBD_stem_30_native_res5,
                      NC_crown_30_native_res5, NC_stem_30_native_res5)
res5_native_30$sample_size<-30

res5_all_samples_native<-rbind(res5_native_20,res5_native_25,res5_native_30)
res5_all_samples_native$status<-"native_species_only"


#Native + nonnative pgls @ res 5


FBD_crown_20_nonnative_res5<-readRDS("pgls_results/FBD_crown_resampled_20_nonnative_concensus_res5.rds")
FBD_crown_20_nonnative_res5<-FBD_crown_20_nonnative_res5[-14]
FBD_crown_20_nonnative_res5$Tree<-"FBD_crown"
FBD_stem_20_nonnative_res5<-readRDS("pgls_results/FBD_stem_resampled_20_nonnative_concensus_res5.rds")
FBD_stem_20_nonnative_res5<-FBD_stem_20_nonnative_res5[-14]
FBD_stem_20_nonnative_res5$Tree<-"FBD_stem"
NC_crown_20_nonnative_res5<-readRDS("pgls_results/NC_crown_resampled_20_nonnative_concensus_res5.rds")
NC_crown_20_nonnative_res5<-NC_crown_20_nonnative_res5[-14]
NC_crown_20_nonnative_res5$Tree<-"NC_crown"
NC_stem_20_nonnative_res5<-readRDS("pgls_results/NC_stem_resampled_20_nonnative_concensus_res5.rds")
NC_stem_20_nonnative_res5<-NC_stem_20_nonnative_res5[-14]
NC_stem_20_nonnative_res5$Tree<-"NC_stem"
res5_nonnative_20<-rbind(FBD_crown_20_nonnative_res5, FBD_stem_20_nonnative_res5,
                         NC_crown_20_nonnative_res5, NC_stem_20_nonnative_res5)
res5_nonnative_20$sample_size<-20


FBD_crown_25_nonnative_res5<-readRDS("pgls_results/FBD_crown_resampled_25_nonnative_concensus_res5.rds")
FBD_crown_25_nonnative_res5<-FBD_crown_25_nonnative_res5[-14]
FBD_crown_25_nonnative_res5$Tree<-"FBD_crown"
FBD_stem_25_nonnative_res5<-readRDS("pgls_results/FBD_stem_resampled_25_nonnative_concensus_res5.rds")
FBD_stem_25_nonnative_res5<-FBD_stem_25_nonnative_res5[-14]
FBD_stem_25_nonnative_res5$Tree<-"FBD_stem"
NC_crown_25_nonnative_res5<-readRDS("pgls_results/NC_crown_resampled_25_nonnative_concensus_res5.rds")
NC_crown_25_nonnative_res5<-NC_crown_25_nonnative_res5[-14]
NC_crown_25_nonnative_res5$Tree<-"NC_crown"
NC_stem_25_nonnative_res5<-readRDS("pgls_results/NC_stem_resampled_25_nonnative_concensus_res5.rds")
NC_stem_25_nonnative_res5<-NC_stem_25_nonnative_res5[-14]
NC_stem_25_nonnative_res5$Tree<-"NC_stem"

res5_nonnative_25<-rbind(FBD_crown_25_nonnative_res5, FBD_stem_25_nonnative_res5,
                         NC_crown_25_nonnative_res5, NC_stem_25_nonnative_res5)
res5_nonnative_25$sample_size<-25



FBD_crown_30_nonnative_res5<-readRDS("pgls_results/FBD_crown_resampled_30_nonnative_concensus_res5.rds")
FBD_crown_30_nonnative_res5<-FBD_crown_30_nonnative_res5[-14]
FBD_crown_30_nonnative_res5$Tree<-"FBD_crown"
FBD_stem_30_nonnative_res5<-readRDS("pgls_results/FBD_stem_resampled_30_nonnative_concensus_res5.rds")
FBD_stem_30_nonnative_res5<-FBD_stem_30_nonnative_res5[-14]
FBD_stem_30_nonnative_res5$Tree<-"FBD_stem"
NC_crown_30_nonnative_res5<-readRDS("pgls_results/NC_crown_resampled_30_nonnative_concensus_res5.rds")
NC_crown_30_nonnative_res5<-NC_crown_30_nonnative_res5[-14]
NC_crown_30_nonnative_res5$Tree<-"NC_crown"
NC_stem_30_nonnative_res5<-readRDS("pgls_results/NC_stem_resampled_30_nonnative_concensus_res5.rds")
NC_stem_30_nonnative_res5<-NC_stem_30_nonnative_res5[-14]
NC_stem_30_nonnative_res5$Tree<-"NC_stem"

res5_nonnative_30<-rbind(FBD_crown_30_nonnative_res5, FBD_stem_30_nonnative_res5,
                         NC_crown_30_nonnative_res5, NC_stem_30_nonnative_res5)
res5_nonnative_30$sample_size<-30

res5_all_samples_nonnative<-rbind(res5_nonnative_20,res5_nonnative_25,res5_nonnative_30)
res5_all_samples_nonnative$status<-"native_nonnative_species"

res5_all<-rbind(res5_all_samples_nonnative, res5_all_samples_native)

res5_all<-res5_all %>% filter(!r2 == 0)

res5_all %>%
  mutate(slope_QWD = as.numeric(slope_QWD),
         sample_size = as.factor(sample_size)) %>%
  ggplot(.) +
  geom_boxplot(aes(x = sample_size, y = slope_QWD, fill = Tree)) +
  facet_wrap(~status)+
  labs(title = "@ 250km2")










