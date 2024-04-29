# predict worker biomass with traits

#read in biomass data

biomass<- read.csv("traits/biomass.csv", encoding = "UTF-8")

biomass_true<-biomass %>% drop_na()
biomass_need<-biomass %>% filter(is.na(biomass_indiv))

#read in worker traits

#Load worker data
worker<-read.csv("traits/worker_traits.csv", encoding = "UTF-8")
#get genus column
worker$genus<-str_extract(worker$SPECIES, "[^\\s]+")


#get genus means
genus_means_w<-worker %>%
  group_by(genus) %>%
  summarise( HW = mean(HW, na.rm =T),
             HL = mean(HL, na.rm =T),
             ML = mean(ML, na.rm =T),
             EL = mean(EL, na.rm =T),
             WL = mean(WL, na.rm =T))
###
worker_traits_sp<-worker %>%
  drop_na(HW) %>%
  group_by(SPECIES) %>%
  summarise(
    HW = mean(HW, na.rm =T),
    HL = mean(HL, na.rm =T),
    ML = mean(ML, na.rm =T),
    EL = mean(EL, na.rm =T),
    WL_min = min(WL, na.rm =T)
  )


#left join

left_join(biomass_true, worker_traits_sp, by = "SPECIES") %>%
  ggplot(.) +
  geom_point(aes(x = WL_min, y = biomass_indiv)) +
  geom_smooth(method = "loess", aes(x = WL_min, y = biomass_indiv)) +
  scale_y_log10()




combined<-left_join(biomass_true, worker_traits_sp, by = "SPECIES")
m1<-lm(log(biomass_indiv)~ poly(WL_min,2), data = combined)
summary(m1)


#leftjoin biomass needed
newdat<-left_join(biomass_need, worker_traits_sp, by = "SPECIES") %>%
  dplyr::select(WL_min)


#now we use the predicted biomass values
biomass_need$biomass_indiv<-exp(predict.lm(m1, newdata = newdat))


new_biomass<-(rbind(biomass_need, biomass_true))
saveRDS(new_biomass,"traits/imputedbiomass.rds")
