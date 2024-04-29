#Supplementary figures 1
#Checking congruence of PGLS across Tree types

#load tree data
#native
FBD_crown_res4_native<-readRDS("pgls_results/FBD_crown_resampled_native_concensus_res4_sensitivity_data.rds")
FBD_crown_res4_native<-FBD_crown_res4_native %>%
  mutate(tree_type = "FBD_crown",
         data_type = "native_species_only")
FBD_stem_res4_native<-readRDS("pgls_results/FBD_stem_resampled_native_concensus_res4_sensitivity_data.rds")
FBD_stem_res4_native<-FBD_stem_res4_native %>%
  mutate(tree_type = "FBD_stem",
         data_type = "native_species_only")
NC_crown_res4_native<-readRDS("pgls_results/NC_crown_resampled_native_concensus_res4_sensitivity_data.rds")
NC_crown_res4_native<-NC_crown_res4_native %>%
  mutate(tree_type = "NC_stem",
         data_type = "native_species_only")

#nonnative
FBD_crown_res4_nonnative<-readRDS("pgls_results/FBD_crown_resampled_nonnative_concensus_res4_sensitivity_data.rds")
FBD_crown_res4_nonnative<-FBD_crown_res4_nonnative %>%
  mutate(tree_type = "FBD_crown",
         data_type = "native_species_nonnative_species")
FBD_stem_res4_nonnative<-readRDS("pgls_results/FBD_stem_resampled_nonnative_concensus_res4_sensitivity_data.rds")
FBD_stem_res4_nonnative<-FBD_stem_res4_nonnative %>%
  mutate(tree_type = "FBD_stem",
         data_type = "native_species_nonnative_species")
NC_crown_res4_nonnative<-readRDS("pgls_results/NC_crown_resampled_nonnative_concensus_res4_sensitivity_data.rds")
NC_crown_res4_nonnative<-NC_crown_res4_nonnative %>%
  mutate(tree_type = "NC_stem",
         data_type = "native_species_nonnative_species")

res_4_all<-rbind(FBD_crown_res4_native,
      FBD_stem_res4_native,
      NC_crown_res4_native)

res_4_all_both<-rbind(FBD_crown_res4_nonnative,
                      FBD_stem_res4_nonnative,
                      NC_crown_res4_nonnative)

ggplot(data= res_4_all) +
  geom_boxplot(aes(x = as.factor(sample_size), y = mean_slope_QWD)) +
  facet_grid(tree_type~data_type)

ggplot(data= res_4_all_both) +
  geom_boxplot(aes(x = as.factor(sample_size), y = mean_slope_QWD)) +
  facet_grid(tree_type~data_type)
