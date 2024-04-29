library(tidyverse)
library(patchwork)

#Figure 3 in manuscript
differences4<-readRDS("difference_data/differences_res4.rds")
differences3<-readRDS("difference_data/differences_res3.rds")
differences5<-readRDS("difference_data/differences_res5.rds")

colnames(differences4[[1]])
colnames(differences3[[1]])
ncol(differences5[[1]])

differences4=lapply(differences4, function(x){
  x<-x%>% mutate(resolution = "~1,770km2")  
  return(x)
})
differences5=lapply(differences5, function(x){
  x<-x%>% mutate(resolution = "~250km2")  
  return(x)
})
differences3=lapply(differences3, function(x){
  x<-x%>% mutate(resolution = "~12,330km2")  
  return(x)
})

Sample_size<-c(20,25,30)
differences4<-map2(.x = differences4, .y = Sample_size, function(x,y){
  x<-x %>%
    mutate(PGLS_sample_size = y)
  return(x)
})

differences5<-map2(.x = differences5, .y = Sample_size, function(x,y){
  x<-x %>%
    mutate(PGLS_sample_size = y)
  return(x)
})

differences3<-map2(.x = differences3, .y = Sample_size, function(x,y){
  x<-x %>%
    mutate(PGLS_sample_size = y)
  return(x)
})

differences4<-do.call(rbind, differences4)
differences3<-do.call(rbind, differences3)
differences5<-do.call(rbind, differences5)

differences_all<-rbind(differences5,differences3, differences4)
differences_all$resolution<-as.factor(differences_all$resolution)
differences_all$PGLS_sample_size<-as.factor(differences_all$PGLS_sample_size)

levels(differences_all$resolution)

levels(differences_all$PGLS_sample_size)
str(differences_all)

differences_all<-differences_all %>%
  mutate(name = fct_relevel(resolution, 
                            "~250 km2", "~1,770 km2", "~12,330 km2"))


levels(differences_all$name)<-c(expression("250"~km^2),
                                expression("1,770"~km^2),
                                expression("12,330"~km^2))
exotic_plot_slope<-ggplot(data = differences_all ) +
  geom_point(aes(x = exotic_perc_avg*100, y= slope_diff,
                 fill = PGLS_sample_size), pch = 21, 
             size = 3,
             alpha = 0.6) +
  geom_smooth(aes(x = exotic_perc_avg*100, y= slope_diff,
                  color = PGLS_sample_size,
                  fill = PGLS_sample_size),alpha = 0.2, method = "lm") +
  geom_hline(yintercept = 0, lty =2) +
  facet_wrap(~name, labeller = label_parsed) +
  theme_light() +
  labs(x = "Nonnative ant %",
       color = "PGLS sample size",
       y = "Slope Difference", fill = "PGLS sample size") +
  scale_fill_viridis_d(option = "F") +
  scale_color_viridis_d(option = "F") +
  theme(axis.title.x = element_blank())
exotic_plot_slope


exotic_plot_r2<-ggplot(data = differences_all ) +
  geom_point(aes(x = exotic_perc_avg*100, y= r2_diff,
                 fill = PGLS_sample_size), pch = 21, color = "black", 
             size = 3,
             alpha = 0.6) +
  geom_smooth(aes(x = exotic_perc_avg*100, y= r2_diff,
                  color = PGLS_sample_size,
                  fill = PGLS_sample_size),alpha = 0.2, method = "lm") +
  geom_hline(yintercept = 0, lty =2) +
  
  facet_wrap(~name,labeller = label_parsed) +
  theme_light() +
  labs(x = "Nonnative ant %",
       color = "PGLS sample size",
       y = "R2 Difference", fill = "PGLS sample size") +
  scale_fill_viridis_d(option = "F")+
  scale_color_viridis_d(option = "F") +
  theme(axis.title.x = element_blank())


exotic_int<-ggplot(data = differences_all ) +
  geom_point(aes(x = exotic_perc_avg*100, y= int_diff,
                 fill = PGLS_sample_size), pch = 21, color = "black", 
             size = 3,
             alpha = 0.6) +
  geom_smooth(aes(x = exotic_perc_avg*100, y= int_diff,
                  color = PGLS_sample_size,
                  fill = PGLS_sample_size),alpha = 0.2, method = "lm") +
  geom_hline(yintercept = 0, lty =2) +
  
  facet_wrap(~name,labeller = label_parsed) +
  theme_light() +
  labs(x = "Nonnative ant %",
       color = "PGLS sample size",
       y = "Intercept Difference", fill = "PGLS sample size") +
  scale_fill_viridis_d(option = "F") +
  scale_color_viridis_d(option = "F") 

exotic_plot_slope + exotic_plot_r2 +exotic_int +
  plot_layout(nrow = 3, guides = "collect", axis_titles = "collect")




#Figure 4
#now plot model estimates
model_est=read.csv("manuscript_data_model_plot.csv")

model_est<-model_est %>%
  mutate(Grain = fct_relevel(Grain, 
                            "~250km2", "~1,770km2", "~12,330km2"))%>%
  mutate(Predictor = fct_relevel(Predictor, 
                             "ATR", "DLC", "ATR:DLC")) 



levels(model_est$Grain)<-c(expression("250"~km^2),
                                expression("1,770"~km^2),
                                expression("12,330"~km^2))
ggplot() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_linerange(data = model_est, 
               aes(x = Effect,xmin = Effect-Effect_sd, xmax = Effect+Effect_sd,
                   y = Predictor, 
                   color = as.factor(sample_size)),
               size = 1,
               position = position_dodge2(width  = 0.8)) +
  geom_point(data = model_est , 
             aes(x = Effect,
                 y = Predictor, 
                 fill = as.factor(sample_size)),
             size = 3,pch =21,
             position = position_dodge2(width = 0.8)) +

  facet_wrap(~Grain,
             nrow = 3,label = label_parsed) +
  theme_dark() +
  labs(x = "Model coefficient", color = "PGLS\nsample\nsize",
       fill = "PGLS\nsample\nsize") +
  theme(axis.title.y = element_blank()) +
  scale_fill_viridis_d(option = "F") +
  scale_color_viridis_d(option = "F")
