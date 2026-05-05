#this script reruns part of the PGLS analysis to gain different metrics based on reviewer requests


#resampling to account for sp richness
#here we rarefy data across many iterations to 
#remove species richness as a potential bias to 
#the slope values
library(sf)
library(ape)
library(nlme)
library(h3jsr)
library(tidyverse)
library(exactextractr)
#load species data
plot<-readRDS("traits/species_final_traits.rds")

remove_species = plot %>%
  filter(QWD_n<0 | is.na(QWD_n)) 

remove_species$SPECIES2<-gsub(" ", ".", remove_species$SPECIES)

plot<-plot %>% filter(!SPECIES %in% remove_species$SPECIES) %>%
  drop_na(QWD_n)
plotspecies<-gsub(" ", ".", plot$SPECIES)
plot$Status

#read in point data for the three states

#read in point data for the three states
intersectionsh3<-readRDS("geo_data/points.rds")

intersectionsh3<-intersectionsh3 %>%
  st_as_sf(coords = c("lon", "lat"),
           crs =  "EPSG:4326")

#read in data from Ohyama et al.
intersectionsh4<-readRDS("geo_data/cleaned_ant_point_data.rds")
intersectionsh4<-intersectionsh4 %>%
  st_as_sf(coords = c("LocLongitude", "LocLatitude"), crs = st_crs(intersectionsh3))

intersectionsh3<-intersectionsh3 %>%
  filter(species %in% plotspecies) %>%
  dplyr::select(species) %>%
  st_as_sf()

intersectionsh4<-intersectionsh4 %>%
  filter(Genus_species %in% plotspecies) %>%
  dplyr::select(Genus_species) %>%
  rename(species = Genus_species)

intersectionsh3 = rbind(intersectionsh3,intersectionsh4)

data<-readRDS("geo_data/GBIF_clean.rds")

intersection5<-data %>%
  select(species, decimalLatitude, decimalLongitude) %>%
  drop_na() %>%
  st_as_sf(coords = c("decimalLongitude","decimalLatitude"),
           crs= st_crs(intersectionsh3))

intersectionsh3 = rbind(intersectionsh3,intersection5)




#using uber h3 to get h3 addresses for resoltution 4
all_res<-point_to_cell(intersectionsh3, res = 4, simple = FALSE)
#get native nonnative percentages 

left_join(all_res, plot %>% select(SPECIES, Status), by =c("species" ="SPECIES")) %>%
  group_by(res, Status) %>%
  summarise(abund = n(),
            n_sp = n_distinct(species)) %>%
  ungroup() %>%
  pivot_wider(id_cols = res, names_from = Status, values_from = c(abund, n_sp)) %>%
  rowwise() %>%
  mutate(total_sp = sum(n_sp_n, n_sp_e, na.rm = T))

left_join(all_res, plot %>% select(SPECIES, Status), by =c("species" ="SPECIES")) %>%
  group_by(res) %>%
  summarise(abund = n(),
            n_sp = n_distinct(species))

#quick code here to rename colnames and fix genus names 
colnum<-which(str_detect(colnames(all_res), "resolution"))
colnames(all_res)[colnum]<-"res"
all_res$species<-gsub("\\."," ", all_res$species)


#filter for only species we have in plot
all_res<-all_res %>%
  filter(species %in% plot$SPECIES)


#get a list of the species we used imputed data for
imputed_sp<-plot %>%
  filter(QWDimpute == "IMPUTED")

topimputedgenera<-imputed_sp %>%
  mutate(genus = str_extract(SPECIES, "[^\\s]+")) %>%
  group_by(genus) %>%
  summarise()
  
total_comp<-all_res %>%
  group_by(species) %>%
  count() %>%
  mutate(imputed = ifelse(species %in% imputed_sp$SPECIES, "imputed",NA)) %>%
  ungroup() %>%
  mutate(perc = 100 *( n/sum(n)))

total_comp %>%
  group_by(imputed) %>%
  summarise(perc= sum(perc))

total_comp<-all_res %>%
  group_by(res,species) %>%
  count() %>%
  mutate(imputed = ifelse(species %in% imputed_sp$SPECIES, "imputed",NA)) %>%
  group_by(res) %>%
  mutate(perc = 100 *( n/sum(n)))

total_n = total_comp %>%
  group_by(res) %>%
  summarise(n_sp = n(),
            abundance = sum(n)) %>%
  ungroup()

#calculate per spatial unit of imputeted species %
perc_by_Res<-total_comp %>%
  group_by(res,imputed) %>%
  summarise(perc= sum(perc)) %>% 
  filter(imputed  == "imputed") %>%
  left_join(.,total_n, by = "res") %>%
  filter(n_sp > 15)


sinvicta<-total_comp %>% 
  filter(species == "Solenopsis invicta")

saveRDS()
sinvicta %>%
  group_by(res) %>%
  summarise(perc= sum(perc))
#now get S for each res and filter for cells with greater than 15 species
#remove exotic species
S_cells_20<-all_res %>% 
  group_by(res) %>%
  filter(!duplicated(species)) %>%
  summarise(S = n()) %>%
  ungroup() %>%
  filter(S>19)



S_cells_25<-all_res %>% 
  group_by(res) %>%
  filter(!duplicated(species)) %>%
  summarise(S = n()) %>%
  ungroup() %>%
  filter(S>24)


S_cells_30<-all_res %>% 
  group_by(res) %>%
  filter(!duplicated(species)) %>%
  summarise(S = n()) %>%
  ungroup() %>%
  filter(S>29)


#now get a community roster for every cell

rosters_20<-all_res %>% 
  filter(res %in% S_cells_20$res) %>%
  group_by(res) %>%
  filter(!duplicated(species)) %>%
  nest()



rosters_20 %>% unnest() %>%
  left_join(., plot, by = c("species" = "SPECIES")) %>%
  group_by(res) %>%
  summarise(S = n(),
            QWD = mean(QWD_n)) %>%
  ggplot(.) +
  geom_point(aes(x = QWD, y = S))



rosters_25<-all_res %>% 
  filter(res %in% S_cells_25$res) %>%
  group_by(res) %>%
  filter(!duplicated(species)) %>%
  nest()

rosters_30<-all_res %>% 
  filter(res %in% S_cells_30$res) %>%
  group_by(res) %>%
  filter(!duplicated(species)) %>%
  nest()
#now we resample 20,25,30 species at random, without replacement
# for every cell

rosters_all_20<-list()

for(i in 1:100){
  rosters2<-lapply(seq_along(rosters_20$data), function(x){
    sp_selected<-sample(1:nrow(rosters_20$data[[x]]), 20)
    sp = rosters_20$data[[x]]$species[sp_selected]
    roster<-data.frame(species = sp,
                       res = rosters_20$res[x])
    
    return(roster)
    
  })
  rosters_all_20[[i]]<-do.call(rbind, rosters2)
  
}

#join with native_nonnative to get percentages from simulated data
wtf20<-lapply(rosters_all_20, function(x){
  df<-left_join(x, plot %>% select(SPECIES, Status), by =c("species" ="SPECIES")) %>%
    group_by(res, Status) %>%
    summarise(abund = n(),
              n_sp = n_distinct(species)) %>%
    ungroup() %>%
    pivot_wider(id_cols = res, names_from = Status, values_from = c(abund, n_sp)) %>%
    rowwise() %>%
    mutate(total_sp = sum(n_sp_n, n_sp_e, na.rm = T),
           exotic_perc = n_sp_e/total_sp)
  
  average_exotic = mean(df$exotic_perc, na.rm = T)
  min = min(df$exotic_perc, na.rm = T)
  max = max(df$exotic_perc, na.rm = T)
  
  data.frame(avg_exotic = average_exotic,min = min, max = max)
})

#averages and min and max
wtf20v2df<-do.call(rbind,wtf20)
mean(wtf20v2df$avg_exotic)
mean(wtf20v2df$min)
mean(wtf20v2df$max)

rosters_all_20[[1]] %>% 
  left_join(., plot, by = c("species" = "SPECIES")) %>%
  group_by(res) %>%
  summarise(S = n(),
            QWD = mean(QWD_n)) %>%
  ggplot(.) +
  geom_point(aes(x = QWD, y = S))





rosters_all_25<-list()

for(i in 1:100){
  rosters2<-lapply(seq_along(rosters_25$data), function(x){
    sp_selected<-sample(1:nrow(rosters_25$data[[x]]), 25)
    sp = rosters_25$data[[x]]$species[sp_selected]
    roster<-data.frame(species = sp,
                       res = rosters_25$res[x])
    
    return(roster)
    
  })
  rosters_all_25[[i]]<-do.call(rbind, rosters2)
  
}

wtf25<-lapply(rosters_all_25, function(x){
  df<-left_join(x, plot %>% select(SPECIES, Status), by =c("species" ="SPECIES")) %>%
    group_by(res, Status) %>%
    summarise(abund = n(),
              n_sp = n_distinct(species)) %>%
    ungroup() %>%
    pivot_wider(id_cols = res, names_from = Status, values_from = c(abund, n_sp)) %>%
    rowwise() %>%
    mutate(total_sp = sum(n_sp_n, n_sp_e, na.rm = T),
           exotic_perc = n_sp_e/total_sp)
  
  average_exotic = mean(df$exotic_perc, na.rm = T)
  min = min(df$exotic_perc, na.rm = T)
  max = max(df$exotic_perc, na.rm = T)
  
  data.frame(avg_exotic = average_exotic,min = min, max = max)
})

#averages and min and max
wtf25v2df<-do.call(rbind,wtf25)
mean(wtf25v2df$avg_exotic)
mean(wtf25v2df$min)
mean(wtf25v2df$max)


rosters_all_30<-list()

for(i in 1:100){
  
  rosters2<-lapply(seq_along(rosters_30$data), function(x){
    sp_selected<-sample(1:nrow(rosters_30$data[[x]]), 30)
    sp = rosters_30$data[[x]]$species[sp_selected]
    roster<-data.frame(species = sp,
                       res = rosters_30$res[x])
    
    return(roster)
    
  })
  rosters_all_30[[i]]<-do.call(rbind, rosters2)
}

wtf30<-lapply(rosters_all_30, function(x){
  df<-left_join(x, plot %>% select(SPECIES, Status), by =c("species" ="SPECIES")) %>%
    group_by(res, Status) %>%
    summarise(abund = n(),
              n_sp = n_distinct(species)) %>%
    ungroup() %>%
    pivot_wider(id_cols = res, names_from = Status, values_from = c(abund, n_sp)) %>%
    rowwise() %>%
    mutate(total_sp = sum(n_sp_n, n_sp_e, na.rm = T),
           exotic_perc = n_sp_e/total_sp)
  
  average_exotic = mean(df$exotic_perc, na.rm = T)
  min = min(df$exotic_perc, na.rm = T)
  max = max(df$exotic_perc, na.rm = T)
  
  data.frame(avg_exotic = average_exotic,min = min, max = max)
})

#averages and min and max
wtf30v2df<-do.call(rbind,wtf30)
mean(wtf30v2df$avg_exotic)
mean(wtf30v2df$min)
mean(wtf30v2df$max)


pgls_sets_20<-lapply(rosters_all_20, function(x){
  pgls_data<-x %>%
    left_join(., plot, by =c("species" = "SPECIES")) %>%
    group_by(res) %>%
    nest() 
  pgls_data$ID = 1:nrow(pgls_data)
  return(pgls_data)
})

pgls_data_20<-rosters_all_20[[1]] %>%
  left_join(., plot, by =c("species" = "SPECIES")) %>%
  group_by(res) %>%
  nest() 



pgls_sets_25<-lapply(rosters_all_25, function(x){
  pgls_data<-x %>%
    left_join(., plot, by =c("species" = "SPECIES")) %>%
    group_by(res) %>%
    nest() 
  pgls_data$ID = 1:nrow(pgls_data)
  return(pgls_data)
})

pgls_data_25<-rosters_all_25[[1]] %>%
  left_join(., plot, by =c("species" = "SPECIES")) %>%
  group_by(res) %>%
  nest() 


pgls_sets_30<-lapply(rosters_all_30, function(x){
  pgls_data<-x %>%
    left_join(., plot, by =c("species" = "SPECIES")) %>%
    group_by(res) %>%
    nest() 
  pgls_data$ID = 1:nrow(pgls_data)
  return(pgls_data)
})

pgls_data_30<-rosters_all_30[[1]] %>%
  left_join(., plot, by =c("species" = "SPECIES")) %>%
  group_by(res) %>%
  nest() 



pgls_data_20$ID<-1:nrow(pgls_data_20)
pgls_data_25$ID<-1:nrow(pgls_data_25)
pgls_data_30$ID<-1:nrow(pgls_data_30)
