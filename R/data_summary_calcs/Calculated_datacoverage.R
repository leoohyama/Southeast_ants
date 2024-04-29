#check geographic coverage
#load species data
plot<-readRDS("traits/species_final_traits.rds")

remove_species = plot %>%
  filter(QWD_n<0 | is.na(QWD_n)) 

remove_species$SPECIES2<-gsub(" ", ".", remove_species$SPECIES)

plot<-plot %>% filter(!SPECIES %in% remove_species$SPECIES) %>%
  drop_na(QWD_n)

#load total SE ant list
SE_ants<-read.csv("all_se_ants_final.csv")

#check what percent of list is completed with trait data

not_included<-SE_ants %>%
  filter(!SPECIES %in% plot$SPECIES)




#read in point data for the three states to assess geogrpahic data coverage

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
  dplyr::select(species) %>%
  st_as_sf()

intersectionsh4<-intersectionsh4 %>%
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

#quick code here to rename colnames and fix genus names 
colnum<-which(str_detect(colnames(all_res), "resolution"))
colnames(all_res)[colnum]<-"res"
all_res$species<-gsub("\\."," ", all_res$species)

total<-all_res %>%
  filter(species %in% SE_ants$SPECIES)
what_we_use<-total %>%
  filter(species %in% plot$SPECIES)

#geospatial data not used
total_not_used<-all_res %>%
  filter(!species %in% SE_ants$SPECIES)

total_not_used %>%
  group_by(species) %>%
  summarise(n = n()) %>%
  arrange(desc(n))
