library(tidyverse)
library(ggrepel)


#read in poly
poly<-read.csv("traits/polymorphism.csv",
               encoding = "UTF-8")
#read in status
status<-read.csv("traits/status.csv",
               encoding = "UTF-8")

#read in colony size data
ants_fl<-read.csv(file = "Colony_size_data/ants_Florida_colony.csv",
                  na.strings=c("NA","NaN", " ", ""))
ants_fl<-ants_fl %>%
  mutate(state = "fl") 

ants_al<-read.csv(file = "Colony_size_data/ants_alabama_colony.csv",
                  na.strings=c("NA","NaN", " ", ""))
ants_al<-ants_al %>%
  mutate(state = "al")

ants_ga<-read.csv(file = "Colony_size_data/ants_Georgia_colony.csv",
                  na.strings=c("NA","NaN", " ", ""))
ants_ga<-ants_ga %>%
  mutate(state = "ga")

colony_size<-rbind(ants_al, ants_fl, ants_ga)

#edit names to match
edit<-which(colony_size$SPECIES == "Camponotus pennsylvanicus ")
colony_size[edit,]<-"Camponotus pennsylvanicus"
rm(edit)

sizes<-as.data.frame(unique(colony_size[,c("SPECIES", "size")]))
sizes$size<-as.numeric(sizes$size)
sp_size<-sizes %>%
  drop_na(size) %>%
  group_by(SPECIES) %>%
  summarise(max_size = max(size),
            median_size = median(size),
            sample = n())

write.csv(sp_size, "traits/species_colony_size.csv")


ants_se<-read.csv(file = "Colony_size_data/southeast_Rest_ofants.csv",
                  na.strings=c("NA","NaN", " ", ""))

ants_se<-ants_se %>%
  mutate(state = "se") 
all_antsv2<-rbind(ants_al, ants_fl, ants_ga,ants_se)

SE_ants<-data.frame(SPECIES = unique(all_antsv2$SPECIES))



#get trait data 
#Load worker data
worker<-read.csv("traits/worker_traits.csv", encoding = "UTF-8")

workernumb<-worker %>%
  group_by(SPECIES) %>%
  summarise(n = n())
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
worker_traits_sp<-worker %>%group_by(SPECIES) %>%
  summarise(
    HW = mean(HW, na.rm =T),
    HL = mean(HL, na.rm =T),
    ML = mean(ML, na.rm =T),
    EL = mean(EL, na.rm =T),
    WL = mean(WL, na.rm =T)
  )

#plot of hw by colony

HW_colony<-left_join(worker_traits_sp, sp_size, by = "SPECIES") %>%
  drop_na(median_size)
 


#####
#give ,me the number of NAs for each trait
for (i in seq_along(worker_traits_sp)) {
  na<-length(which(is.na(worker_traits_sp[,i])))
  print(c(colnames(worker_traits_sp)[[i]], na))
}
worker_traits_sp$genus<-str_extract(worker_traits_sp$SPECIES, "[^\\s]+")

#read in queen traits
queen<-read.csv("traits/queen_traits.csv", encoding = "UTF-8")
#get genus column
queen$genus<-str_extract(queen$SPECIES, "[^\\s]+")



#get genus means
genus_means_q<-queen %>%
  group_by(genus) %>%
  summarise( HW = mean(HW, na.rm =T),
             HL = mean(HL, na.rm =T),
             ML = mean(ML, na.rm =T),
             EL = mean(EL, na.rm =T),
             WL = mean(WL, na.rm =T))

queen_traits_sp<-queen %>%group_by(SPECIES) %>%
  summarise(
    HW = mean(HW, na.rm =T),
    HL = mean(HL, na.rm =T),
    ML = mean(ML, na.rm =T),
    EL = mean(EL, na.rm =T),
    WL = mean(WL, na.rm =T)
  )

#give ,me the number of NAs for each trait
for (i in seq_along(queen_traits_sp)) {
  na<-length(which(is.na(queen_traits_sp[,i])))
  print(c(colnames(queen_traits_sp)[[i]], na))
}

#get QWD by each worker specimen

#first impute queen sp data for hw and wl
queen_traits_sp2<-queen_traits_sp %>%
  select(SPECIES, HW,WL) %>%
  mutate(genus = str_extract(SPECIES, "[^ ]+"))

genus_means_q2<-genus_means_q %>%
  select(genus, HW,WL)
queen_traits_sp2$qwdimpute<-NA

for (i in 1:nrow(queen_traits_sp2)) {
  if(is.na(queen_traits_sp2$HW[i])){
    queen_traits_sp2$HW[i]<-genus_means_q2$HW[genus_means_q2$genus == queen_traits_sp2$genus[i]]

  }
  if(is.na(queen_traits_sp2$WL[i])){
    queen_traits_sp2$WL[i]<-genus_means_q2$WL[genus_means_q2$genus == queen_traits_sp2$genus[i]]
    queen_traits_sp2$qwdimpute[i]<-"IMPUTED"
  }
  
  
}

variationQWD<-left_join(worker, queen_traits_sp2 %>%
            select(-genus), by = "SPECIES") %>%
  select(SPECIES, HW.x, HW.y, WL.x, WL.y) %>%
  mutate(QWD_n = 100 * ( (2 * (WL.y - WL.x)) / (WL.y + WL.x)),
         QWD_hw_n = 100 * ( (2 * (HW.y - HW.x)) / (HW.y + HW.x))) %>%
  group_by(SPECIES) %>%
  summarise(stdevQWD = sd(QWD_n, na.rm = T),
            stdevQWD_h = sd(QWD_hw_n, na.rm = T))

#get genus QWD_var
variationQWD_genus<-variationQWD %>%
  mutate(genus = str_extract(SPECIES, "[^ ]+")) %>%
  group_by(genus) %>%
  summarise(stdevQWD_n_g = mean(stdevQWD, na.rm =T),
            stdevQWD_h_g = mean(stdevQWD_h, na.rm =T))

#add genus means for nas
variationQWD
variationQWD$qwdsdevimpute<-NA
for( i in 1:nrow(variationQWD)){
  if(is.na(variationQWD$stdevQWD[i])){
    variationQWD$stdevQWD[i] = variationQWD_genus$stdevQWD_n_g[variationQWD_genus$genus == str_extract(variationQWD$SPECIES[i], "[^ ]+")]
    variationQWD$qwdsdevimpute[i]<-"IMPUTED"
  }
  if(is.na(variationQWD$stdevQWD_h[i])){
    variationQWD$stdevQWD_h[i] =variationQWD_genus$stdevQWD_h_g[variationQWD_genus$genus == str_extract(variationQWD$SPECIES[i],
                                                                                                        "[^ ]+")]
  }
}

variationQWD_corrected<-variationQWD

#get genus QWD for imputation 
QWD_traits_impute<-left_join(genus_means_q, genus_means_w, by = "genus")

QWD_traits_impute_genus<-QWD_traits_impute %>%
  mutate(QWD = (WL.x/WL.y)*100,
         QWD_n = 100 * ( (2 * (WL.x - WL.y)) / (WL.x + WL.y)),
         QWD_HL = 100 * ( (2 * (HL.x - HL.y)) / (HL.x + HL.y)),
         QWD_EL = 100 * ( (2 * (EL.x - EL.y)) / (EL.x + EL.y)),
         QWD_ML = 100 * ( (2 * (ML.x - ML.y)) / (ML.x + ML.y)),
         QWD_HW = 100 * ( (2 * (HW.x - HW.y)) / (HW.x + HW.y))
         
  ) %>%
  dplyr::select(genus, QWD,QWD_n, QWD_HW)


#join by queen
QWD_traits2<-left_join(queen_traits_sp, worker_traits_sp, by = "SPECIES")

QWD_traits<-QWD_traits2 %>%
  mutate(QWD = (WL.x/WL.y)*100,
         QWD_n = 100 * ( (2 * (WL.x - WL.y)) / (WL.x + WL.y)),
         QWD_HL = 100 * ( (2 * (HL.x - HL.y)) / (HL.x + HL.y)),
         QWD_EL = 100 * ( (2 * (EL.x - EL.y)) / (EL.x + EL.y)),
         QWD_ML = 100 * ( (2 * (ML.x - ML.y)) / (ML.x + ML.y)),
         QWD_HW = 100 * ( (2 * (HW.x - HW.y)) / (HW.x + HW.y))
         
  ) %>%
  dplyr::select(SPECIES,genus, QWD,QWD_n, QWD_HW)

#now impute QWD using genus means
#impute missing data with means
QWD_traits<-QWD_traits %>% drop_na(genus)
QWD_traits$QWDimpute<-NA
for(i in 1:nrow(QWD_traits)){


  if(is.na(QWD_traits$QWD_HW[[i]])){
    QWD_traits$QWD_HW[[i]] = QWD_traits_impute_genus$QWD_HW[QWD_traits_impute_genus$genus == QWD_traits$genus[[i]]]
  }

  if(is.na(QWD_traits$QWD_n[[i]])){
    QWD_traits$QWD_n[[i]] = QWD_traits_impute_genus$QWD_n[QWD_traits_impute_genus$genus == QWD_traits$genus[[i]]]
    QWD_traits$QWDimpute[[i]]<-"IMPUTED"
  }else{
    QWD_traits$QWDimpute[[i]]<-"NOT_IMPUTED"
  }
  if(is.na(QWD_traits$QWD[[i]])){
    QWD_traits$QWD[[i]] = QWD_traits_impute_genus$QWD[QWD_traits_impute_genus$genus == QWD_traits$genus[[i]]]
  }
}



#now join with colony
QWD_only<-QWD_traits %>% dplyr::select(SPECIES, QWD_HW, QWD, QWD_n)
plot<-left_join(sp_size, QWD_only, by = "SPECIES")
plot<-left_join(plot, variationQWD_corrected, by = "SPECIES")

plot<-left_join(plot, poly, by = "SPECIES")
plot<-left_join(plot, status, by = "SPECIES")
plot<-plot %>% drop_na(c("Status", "poly_id"))

biomass<-readRDS("traits/imputedbiomass.rds")
plot<-left_join(plot, biomass, by = "SPECIES") %>%
  mutate(biomass = biomass_indiv * median_size,
         biomass_max = biomass_indiv * max_size)

saveRDS(plot, "traits/species_final_traits.rds")


QWD_only<-QWD_traits %>% dplyr::select(SPECIES, QWD_HW, QWD, QWD_n, QWDimpute)
plot<-left_join(sp_size, QWD_only, by = "SPECIES")
plot<-left_join(plot, variationQWD_corrected, by = "SPECIES")

plot<-left_join(plot, poly, by = "SPECIES")
plot<-left_join(plot, status, by = "SPECIES")
plot<-plot %>% drop_na(c("Status", "poly_id"))

biomass<-readRDS("traits/imputedbiomass.rds")
plot<-left_join(plot, biomass, by = "SPECIES") %>%
  mutate(biomass = biomass_indiv * median_size,
         biomass_max = biomass_indiv * max_size)

data_final<-left_join(plot, workernumb, by = 'SPECIES') %>%
  filter(!QWD_n<0) %>%
  filter(!is.na(QWD_n))

data_final %>%
  group_by(QWDimpute) %>%
  summarise(n=n())

##nonnative plot
plot<-readRDS( "traits/species_final_traits.rds")

ggplot() +
  geom_point(data = plot %>%
             filter(QWD_n >0), aes(x = QWD_n, y = median_size),
             pch =21, fill = "grey", size = 3) +
  geom_jitter(data = plot %>%
               filter(QWD_n >0) %>%
               filter(SPECIES %in% c(
                                     "Pheidole dentata",
                                     "Monomorium minimum",
                                     "Pheidole morrisii", 
                                     "Solenopsis geminata")), aes(x = QWD_n, y = median_size),
             pch =21, fill = "red", size = 4) +
  geom_point(data = plot %>%
               filter(QWD_n >0) %>%
               filter(SPECIES %in% c(
                 "Solenopsis invicta")), aes(x = QWD_n, y = median_size),
             pch =21, fill = "black", size = 4) +
  geom_text_repel(data = plot %>%
                    filter(QWD_n >0) %>%
                    filter(SPECIES %in% c(
                      "Pheidole dentata",
                      "Monomorium minimum", "Pheidole morrisii",
                      "Solenopsis geminata",
                      "Solenopsis invicta")), 
                  aes(x = QWD_n, y = median_size, label = SPECIES),
                  size = 5,
                  force             = 0.5,
                  nudge_x           = 0,
                  direction         = "y",
                  hjust             = 0.5,
                  vjust = 1.5,
                  segment.size      = 0.3,
                  segment.curvature = 0
  ) +
  scale_y_log10() +
  theme_linedraw() +
  labs(x = "Queen-worker dimorphism (%)", y = "Colony size (log-scale)") +
  theme(axis.title = element_text(size = 16, face = "bold"))



