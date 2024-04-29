#this is a function that get ,model parameters by spatial assemblage for
#sensitivity analysis
library(tidyverse)
#samplesize<-20
#x<-resampled_20
#model_param %>% filter(res == "8444eb7ffffffff")

sense_data<-function(x, samplesize){

  
  res<-data.frame(res = unique(x$res))
  
  null_number<-x %>%
    filter(p_QWD == "null") %>%
    group_by(res) %>%
    summarise(n_null = n()) 
  
  null_number<-left_join(res, null_number, by = "res") %>%
    mutate(sample_size = samplesize) %>%
    mutate(n_null = ifelse(is.na(n_null), 0, n_null))
  
  p_val_pass_number<-x %>%
    filter(!p_QWD == "null") %>%
    mutate(p_QWD = as.numeric(p_QWD)) %>%
    group_by(res) %>%
    filter(p_QWD<0.05) %>%
    summarise(pval_pass = n()) 
  
  updated_tb<-left_join(null_number, p_val_pass_number,by = "res") %>%
    mutate(pval_pass = ifelse(is.na(pval_pass), 0, pval_pass))
  
  normality_checks<-x %>%
    filter(!p_QWD == "null") %>%
    filter(normality>0.05) %>%
    group_by(res) %>%
    summarise(normality_pass = n())
  
  updated_tb2<-left_join(updated_tb, normality_checks,by = "res") %>%
    mutate(normality_pass = ifelse(is.na(normality_pass), 0, normality_pass))
  
  homogeneity_checks<-x %>%
    filter(!p_QWD == "null") %>%
    filter(homogeneity>0.05) %>%
    group_by(res) %>%
    summarise(homogeniety_pass = n())
  
  updated_tb3<-left_join(updated_tb2, homogeneity_checks,by = "res") %>%
    mutate(homogeniety_pass = ifelse(is.na(homogeniety_pass), 0, homogeniety_pass))
  
  both_checks<-x %>%
    filter(!p_QWD == "null") %>%
    filter(homogeneity>0.05 & normality > 0.05) %>%
    group_by(res) %>%
    summarise(both_pass = n())
  
  updated_tb4<-left_join(updated_tb3, both_checks,by = "res") %>%
    mutate(both_pass = ifelse(is.na(both_pass), 0, both_pass))


  if("exotic_perc" %in% colnames(resampled_20)){
    model_param<-x %>%
      filter(!p_QWD == "null") %>%
      filter(! r2 == "error") %>%
      group_by(res) %>%
      mutate_at(vars(-2), as.numeric) %>%
      filter(normality > 0.05 & homogeneity > 0.05) %>%
      summarise(
        model_used_no_sig = n(),
        exotic_perc_avg = mean(exotic_perc, na.rm=T),
        med_pQWD = median(p_QWD, na.rm=T),
        med_psdQWD = median(p_sdQWD, na.rm=T),
        mean_slope_QWD = mean(slope_QWD, na.rm=T),
        mean_slope_error_QWD = mean(slope_error_QWD, na.rm =T),
        mean_slope_error_sdQWD = mean(slope_error_sdQWD, na.rm= T),
        mean_slope_sdQWD = mean(slope_sdQWD, na.rm=T),
        mean_intercept = mean(intercept, na.rm=T),
        mean_r2 = mean(r2, na.rm =T)
      ) 
  }else{
  model_param<-x %>%
    filter(!p_QWD == "null") %>%
    filter(! r2 == "error") %>%
    group_by(res) %>%
    mutate_at(vars(-2), as.numeric) %>%
    filter(normality > 0.05 & homogeneity > 0.05) %>%
    summarise(
      model_used_no_sig = n(),
      med_pQWD = median(p_QWD, na.rm=T),
      med_psdQWD = median(p_sdQWD, na.rm=T),
      mean_slope_QWD = mean(slope_QWD, na.rm=T),
      mean_slope_error_QWD = mean(slope_error_QWD, na.rm =T),
      mean_slope_error_sdQWD = mean(slope_error_sdQWD, na.rm= T),
      mean_slope_sdQWD = mean(slope_sdQWD, na.rm=T),
      mean_intercept = mean(intercept, na.rm=T),
      mean_r2 = mean(r2, na.rm =T)
    ) 
  }
  updated_tb5<-left_join(updated_tb4, model_param,by = "res") 
  
  if("exotic_perc" %in% colnames(resampled_20)){
  model_param_sig<-x %>%
    filter(!p_QWD == "null") %>%
    filter(! r2 == "error") %>%
    mutate(p_QWD = as.numeric(p_QWD)) %>%
    filter(p_QWD < 0.05) %>%
    group_by(res) %>%
    mutate_at(vars(-2), as.numeric) %>%
    filter(normality > 0.05 & homogeneity > 0.05) %>%
    summarise(
      model_used_sig = n(),
      exotic_perc_avg_sig = mean(exotic_perc, na.rm=T),
      mean_slope_QWD_sig = mean(slope_QWD, na.rm=T),
      mean_slope_error_QWD_sig = mean(slope_error_QWD, na.rm =T),
      mean_slope_error_sdQWD_sig = mean(slope_error_sdQWD, na.rm= T),
      mean_slope_sdQWD_sig = mean(slope_sdQWD, na.rm=T),
      mean_intercept_sig = mean(intercept, na.rm=T),
      mean_r2_sig = mean(r2, na.rm =T)
    ) }else{
      model_param_sig<-x %>%
        filter(!p_QWD == "null") %>%
        filter(! r2 == "error") %>%
        mutate(p_QWD = as.numeric(p_QWD)) %>%
        filter(p_QWD < 0.05) %>%
        group_by(res) %>%
        mutate_at(vars(-2), as.numeric) %>%
        filter(normality > 0.05 & homogeneity > 0.05) %>%
        summarise(
          model_used_sig = n(),
          mean_slope_QWD_sig = mean(slope_QWD, na.rm=T),
          mean_slope_error_QWD_sig = mean(slope_error_QWD, na.rm =T),
          mean_slope_error_sdQWD_sig = mean(slope_error_sdQWD, na.rm= T),
          mean_slope_sdQWD_sig = mean(slope_sdQWD, na.rm=T),
          mean_intercept_sig = mean(intercept, na.rm=T),
          mean_r2_sig = mean(r2, na.rm =T))
    }
  updated_tb6<-left_join(updated_tb5, model_param_sig,by = "res") 
  
  
 
  return(updated_tb6)
  
}
