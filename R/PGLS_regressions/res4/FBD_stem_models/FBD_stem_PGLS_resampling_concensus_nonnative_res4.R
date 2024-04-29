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
#load trees
FBD_crown<-read.tree("Phylo_data/Pruned_Renamed_Trees/15K_FBD_crown_posterior_LEO_15Nov202.tre")
FBD_stem<-read.tree("Phylo_data/Pruned_Renamed_Trees/15K_FBD_stem_posterior_LEO_15Nov202.tre")
NC_crown<-read.tree("Phylo_data/Pruned_Renamed_Trees/15k_NCuniform_crown_posterior_LEO_15Nov202.tre")
NC_stem<-read.tree("Phylo_data/Pruned_Renamed_Trees/15K_FBD_stem_posterior_LEO_15Nov202.tre")

#we make concensus tree
concensus_tree<-phytools::consensus.edges(FBD_stem)

#load species data
plot<-readRDS("traits/species_final_traits.rds")

remove_species = plot %>%
  filter(QWD_n<0 | is.na(QWD_n)) 

remove_species$SPECIES2<-gsub(" ", ".", remove_species$SPECIES)

plot<-plot %>% filter(!SPECIES %in% remove_species$SPECIES) %>%
  drop_na(QWD_n)
plotspecies<-gsub(" ", ".", plot$SPECIES)
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

#quick code here to rename colnames and fix genus names 
colnum<-which(str_detect(colnames(all_res), "resolution"))
colnames(all_res)[colnum]<-"res"
all_res$species<-gsub("\\."," ", all_res$species)



#filter for only species we have in plot
all_res<-all_res %>%
  filter(species %in% plot$SPECIES)

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



saveRDS(pgls_data_20,"pgls_results/PGLS_datasets_resample_nonnative/pgls_data_res4_s20_FBD_stem.rds")
saveRDS(pgls_data_25,"pgls_results/PGLS_datasets_resample_nonnative/pgls_data_res4_s25_FBD_stem.rds")
saveRDS(pgls_data_30,"pgls_results/PGLS_datasets_resample_nonnative/pgls_data_res4_s30_FBD_stem.rds")

#set the tree
trees<-concensus_tree

#run at 20

results_20<-lapply(seq_along(pgls_sets_20), function(d){
  
  pgls_data<-pgls_sets_20[[d]]
  slopes<-list()
  
  cv<-seq_along(pgls_data$data)
  
  resample_pgls<-lapply(cv, function(x){
    
    df<-as.data.frame(pgls_data$data[[x]])
    df$species<- gsub(" ", ".", df$species)
    word<-(paste(x, pgls_data$res[x], sep ="_"))
    print(paste(d, word, sep = "_"))
    exotic_perc = which(df$Status == "e")
    if(length(exotic_perc) == 0){
      exotic_perc<-0
    }else{
      exotic_perc<-length(which(df$Status == "e"))/nrow(df)
    }
    #set up tree
    tree_pgls<-drop.tip(trees,
                        setdiff(trees$tip.label,
                                df$species))
    
    #run different models
    pagel_1<-tryCatch({gls(log10(median_size) ~ scale(stdevQWD)+ scale(QWD_n), data = df, 
                           correlation=corPagel(1, tree_pgls, form = ~species),
                           method = "ML")}, error = function(err) {return("error")})
    pagel_0_5<-tryCatch({gls(log10(median_size) ~ scale(stdevQWD)+ scale(QWD_n), data = df, 
                             correlation=corPagel(0.5, tree_pgls, form = ~species),
                             method = "ML")}, error = function(err) {return("error")})
    cormartin<-tryCatch({gls(log10(median_size) ~ scale(stdevQWD) + scale(QWD_n), data = df, 
                             correlation=corMartins(1, tree_pgls, form = ~species),
                             method = "ML")}, error = function(err) {return("error")})
    
    Null_model<-tryCatch({lm(log10(median_size) ~ 1,data = df)}, error = function(err) {return("error")})
    
    #make into list
    model_list<-list(Null_model, pagel_0_5,pagel_1, cormartin)
    
    
    #add an if clause to skip if all models fail:
    if(length(which(model_list == "error")) == 4){
      tables<-data.frame(ID = x,
                         res = pgls_data$res[x],
                         resample_no = d,
                         p_QWD = "error",
                         p_sdQWD = "error",
                         r2 = "error",
                         slope_QWD = "error",
                         slope_sdQWD = "error",
                         slope_error_QWD = "error",
                         slope_error_sdQWD = "error",
                         normality = "error",
                         homogeneity = "error",
                         intercept = "error",
                         exotic_perc = "error")
    }else{ 
      #if none of the models had an error
      if(length(which(model_list == "error"))==0){
        aic_ranks<-bbmle::AICctab(Null_model, pagel_0_5,pagel_1, cormartin,weights =T, sort = T,
                                  mnames = 1:4)
        #get top model
        top_model<- attr(aic_ranks, "row.names")[[1]]
        top_mod_num<-str_extract(top_model, "\\d")
        winning_model<-(model_list[[as.numeric(top_mod_num)]])
        
      }else{
        #if there are some model errors, trim those models away
        model_trim = model_list[which(!model_list == "error")]
        aic_ranks<-bbmle::AICctab(model_trim,weights = T, 
                                  sort = T, mnames = c(which(!model_list == "error")))
        print("get top model")
        top_model<- attr(aic_ranks, "row.names")[[1]]
        top_mod_num<-as.numeric(str_extract(top_model, "\\d"))
        winning_model<-(model_list[[(top_mod_num)]])
      }
      
      coeff_names <- names(coef(winning_model))
      
      wtf<-summary(winning_model)
      
      
      r2<-tryCatch({rr2::R2_lik(winning_model)}, error = function(err){
        return("error")})
      #also check residuals of the model
      if(inherits(winning_model, "lm")){
        
        tables<-data.frame(ID = x,
                           res = pgls_data$res[x],
                           resample_no = d, 
                           p_QWD = "null",
                           p_sdQWD = "null",
                           r2 = r2,
                           slope_QWD =  "null",
                           slope_sdQWD = "null",
                           slope_error_QWD = "null",
                           slope_error_sdQWD = "null",
                           normality =  "null",
                           homogeneity =  "null",
                           intercept =   "null",
                           exotic_perc = "null"
        )
      }else{
        #get model coefs
        p_sdQWD = wtf$tTable[2,4]
        p_QWD = wtf$tTable[3,4]
        
        slope_sdQWD = wtf$tTable[2,1]
        slope_error_sdQWD = wtf$tTable[2,2]
        
        slope_QWD= wtf$tTable[3,1]
        slope_error_QWD = wtf$tTable[3,2]
        
        
        intercept =  wtf$tTable[1,1]
        r2 = r2
        
        
        test<- nortest::lillie.test(resid(winning_model, type = "normalized"))
        lilli_test_vector<-test[2][[1]]
        #now test homogenity of variance
        cortests<-cor.test((as.vector(resid(winning_model, type = "normalized"))),
                           as.vector(winning_model$fitted),method = "pearson")
        corfit_vector1<-cortests$p.value
        print("#assign test")
        tables<-data.frame(ID = x,
                           res = pgls_data$res[x],
                           resample_no = d, 
                           p_QWD = p_QWD,
                           p_sdQWD = p_sdQWD,
                           r2 = r2,
                           slope_QWD = slope_QWD, 
                           slope_sdQWD =slope_sdQWD, 
                           slope_error_QWD = slope_error_QWD,
                           slope_error_sdQWD =slope_error_sdQWD,
                           normality = lilli_test_vector,
                           homogeneity = corfit_vector1,
                           intercept =  intercept,
                           exotic_perc =exotic_perc
        )
      }
      
      
      
      
    }
    return(tables)
  })
  
  return(do.call(rbind,resample_pgls))
  
})

resampled_20<-do.call(rbind,results_20)
rownames(resampled_20)<-NULL


saveRDS(resampled_20, "pgls_results/FBD_stem_resampled_20_nonnative_concensus_res4.rds")

#run at 25
results_25<-lapply(seq_along(pgls_sets_25), function(d){
  
  pgls_data<-pgls_sets_25[[d]]
  slopes<-list()
  
  cv<-seq_along(pgls_data$data)
  
  resample_pgls<-lapply(cv, function(x){
    
    
    df<-as.data.frame(pgls_data$data[[x]])
    df$species<- gsub(" ", ".", df$species)
    word<-(paste(x, pgls_data$res[x], sep ="_"))
    print(paste(d, word, sep = "_"))
    exotic_perc = which(df$Status == "e")
    if(length(exotic_perc) == 0){
      exotic_perc<-0
    }else{
      exotic_perc<-length(which(df$Status == "e"))/nrow(df)
    }
    #set up tree
    tree_pgls<-drop.tip(trees,
                        setdiff(trees$tip.label,
                                df$species))
    
    #run different models
    pagel_1<-tryCatch({gls(log10(median_size) ~ scale(stdevQWD)+ scale(QWD_n), data = df, 
                           correlation=corPagel(1, tree_pgls, form = ~species),
                           method = "ML")}, error = function(err) {return("error")})
    pagel_0_5<-tryCatch({gls(log10(median_size) ~ scale(stdevQWD)+ scale(QWD_n), data = df, 
                             correlation=corPagel(0.5, tree_pgls, form = ~species),
                             method = "ML")}, error = function(err) {return("error")})
    cormartin<-tryCatch({gls(log10(median_size) ~ scale(stdevQWD)+ scale(QWD_n), data = df, 
                             correlation=corMartins(1, tree_pgls, form = ~species),
                             method = "ML")}, error = function(err) {return("error")})
    
    Null_model<-tryCatch({lm(log10(median_size) ~ 1,data = df)}, error = function(err) {return("error")})
    
    #make into list
    model_list<-list(Null_model, pagel_0_5,pagel_1, cormartin)
    
    
    #add an if clause to skip if all models fail:
    if(length(which(model_list == "error")) == 4){
      tables<-data.frame(ID = x,
                         res = pgls_data$res[x],
                         resample_no = d,
                         p_QWD = "error",
                         p_sdQWD = "error",
                         r2 = "error",
                         slope_QWD = "error",
                         slope_sdQWD = "error",
                         slope_error_QWD = "error",
                         slope_error_sdQWD = "error",
                         normality = "error",
                         homogeneity = "error",
                         intercept = "error",
                         exotic_perc = "error")
    }else{ 
      #if none of the models had an error
      if(length(which(model_list == "error"))==0){
        aic_ranks<-bbmle::AICctab(Null_model, pagel_0_5,pagel_1, cormartin,weights =T, sort = T,
                                  mnames = 1:4)
        #get top model
        top_model<- attr(aic_ranks, "row.names")[[1]]
        top_mod_num<-str_extract(top_model, "\\d")
        winning_model<-(model_list[[as.numeric(top_mod_num)]])
        
      }else{
        #if there are some model errors, trim those models away
        model_trim = model_list[which(!model_list == "error")]
        aic_ranks<-bbmle::AICctab(model_trim,weights = T, 
                                  sort = T, mnames = c(which(!model_list == "error")))
        print("get top model")
        top_model<- attr(aic_ranks, "row.names")[[1]]
        top_mod_num<-as.numeric(str_extract(top_model, "\\d"))
        winning_model<-(model_list[[(top_mod_num)]])
      }
      
      coeff_names <- names(coef(winning_model))
      
      wtf<-summary(winning_model)
      
      
      r2<-tryCatch({rr2::R2_lik(winning_model)}, error = function(err){
        return("error")})
      #also check residuals of the model
      if(inherits(winning_model, "lm")){
        
        tables<-data.frame(ID = x,
                           res = pgls_data$res[x],
                           resample_no = d, 
                           p_QWD = "null",
                           p_sdQWD = "null",
                           r2 = r2,
                           slope_QWD =  "null",
                           slope_sdQWD = "null",
                           slope_error_QWD = "null",
                           slope_error_sdQWD = "null",
                           normality =  "null",
                           homogeneity =  "null",
                           intercept =   "null",
                           exotic_perc = "null"
        )
      }else{
        #get model coefs
        p_sdQWD = wtf$tTable[2,4]
        p_QWD = wtf$tTable[3,4]
        
        slope_sdQWD = wtf$tTable[2,1]
        slope_error_sdQWD = wtf$tTable[2,2]
        
        slope_QWD= wtf$tTable[3,1]
        slope_error_QWD = wtf$tTable[3,2]
        
        
        intercept =  wtf$tTable[1,1]
        r2 = r2
        
        
        test<- nortest::lillie.test(resid(winning_model, type = "normalized"))
        lilli_test_vector<-test[2][[1]]
        #now test homogenity of variance
        cortests<-cor.test((as.vector(resid(winning_model, type = "normalized"))),
                           as.vector(winning_model$fitted),method = "pearson")
        corfit_vector1<-cortests$p.value
        print("#assign test")
        tables<-data.frame(ID = x,
                           res = pgls_data$res[x],
                           resample_no = d, 
                           p_QWD = p_QWD,
                           p_sdQWD = p_sdQWD,
                           r2 = r2,
                           slope_QWD = slope_QWD, 
                           slope_sdQWD =slope_sdQWD, 
                           slope_error_QWD = slope_error_QWD,
                           slope_error_sdQWD =slope_error_sdQWD,
                           normality = lilli_test_vector,
                           homogeneity = corfit_vector1,
                           intercept =  intercept,
                           exotic_perc = exotic_perc
        )
      }
      
      
      
      
    }
    return(tables)
  })
  
  return(do.call(rbind,resample_pgls))
  
})
resampled_25<-do.call(rbind,results_25)
rownames(resampled_25)<-NULL



saveRDS(resampled_25, "pgls_results/FBD_stem_resampled_25_nonnative_concensus_res4.rds")

#run at 30
results_30<-lapply(seq_along(pgls_sets_30), function(d){
  
  pgls_data<-pgls_sets_30[[d]]
  slopes<-list()
  
  cv<-seq_along(pgls_data$data)
  
  resample_pgls<-lapply(cv, function(x){
    
    
    df<-as.data.frame(pgls_data$data[[x]])
    df$species<- gsub(" ", ".", df$species)
    word<-(paste(x, pgls_data$res[x], sep ="_"))
    print(paste(d, word, sep = "_"))
    exotic_perc = which(df$Status == "e")
    if(length(exotic_perc) == 0){
      exotic_perc<-0
    }else{
      exotic_perc<-length(which(df$Status == "e"))/nrow(df)
    }
    #set up tree
    tree_pgls<-drop.tip(trees,
                        setdiff(trees$tip.label,
                                df$species))
    
    #run different models
    pagel_1<-tryCatch({gls(log10(median_size) ~ scale(stdevQWD)+ scale(QWD_n), data = df, 
                           correlation=corPagel(1, tree_pgls, form = ~species),
                           method = "ML")}, error = function(err) {return("error")})
    pagel_0_5<-tryCatch({gls(log10(median_size) ~ scale(stdevQWD)+ scale(QWD_n), data = df, 
                             correlation=corPagel(0.5, tree_pgls, form = ~species),
                             method = "ML")}, error = function(err) {return("error")})
    cormartin<-tryCatch({gls(log10(median_size) ~ scale(stdevQWD)+ scale(QWD_n), data = df, 
                             correlation=corMartins(1, tree_pgls, form = ~species),
                             method = "ML")}, error = function(err) {return("error")})
    
    Null_model<-tryCatch({lm(log10(median_size) ~ 1,data = df)}, error = function(err) {return("error")})
    
    #make into list
    model_list<-list(Null_model, pagel_0_5,pagel_1, cormartin)
    
    
    #add an if clause to skip if all models fail:
    if(length(which(model_list == "error")) == 4){
      tables<-data.frame(ID = x,
                         res = pgls_data$res[x],
                         resample_no = d,
                         p_QWD = "error",
                         p_sdQWD = "error",
                         r2 = "error",
                         slope_QWD = "error",
                         slope_sdQWD = "error",
                         slope_error_QWD = "error",
                         slope_error_sdQWD = "error",
                         normality = "error",
                         homogeneity = "error",
                         intercept = "error",
                         exotic_perc = "error")
    }else{ 
      #if none of the models had an error
      if(length(which(model_list == "error"))==0){
        aic_ranks<-bbmle::AICctab(Null_model, pagel_0_5,pagel_1, cormartin,weights =T, sort = T,
                                  mnames = 1:4)
        #get top model
        top_model<- attr(aic_ranks, "row.names")[[1]]
        top_mod_num<-str_extract(top_model, "\\d")
        winning_model<-(model_list[[as.numeric(top_mod_num)]])
        
      }else{
        #if there are some model errors, trim those models away
        model_trim = model_list[which(!model_list == "error")]
        aic_ranks<-bbmle::AICctab(model_trim,weights = T, 
                                  sort = T, mnames = c(which(!model_list == "error")))
        print("get top model")
        top_model<- attr(aic_ranks, "row.names")[[1]]
        top_mod_num<-as.numeric(str_extract(top_model, "\\d"))
        winning_model<-(model_list[[(top_mod_num)]])
      }
      
      coeff_names <- names(coef(winning_model))
      
      wtf<-summary(winning_model)
      
      
      r2<-tryCatch({rr2::R2_lik(winning_model)}, error = function(err){
        return("error")})
      #also check residuals of the model
      if(inherits(winning_model, "lm")){
        
        tables<-data.frame(ID = x,
                           res = pgls_data$res[x],
                           resample_no = d, 
                           p_QWD = "null",
                           p_sdQWD = "null",
                           r2 = r2,
                           slope_QWD =  "null",
                           slope_sdQWD = "null",
                           slope_error_QWD = "null",
                           slope_error_sdQWD = "null",
                           normality =  "null",
                           homogeneity =  "null",
                           intercept =   "null",
                           exotic_perc = "null"
        )
      }else{
        #get model coefs
        p_sdQWD = wtf$tTable[2,4]
        p_QWD = wtf$tTable[3,4]
        
        slope_sdQWD = wtf$tTable[2,1]
        slope_error_sdQWD = wtf$tTable[2,2]
        
        slope_QWD= wtf$tTable[3,1]
        slope_error_QWD = wtf$tTable[3,2]
        
        
        intercept =  wtf$tTable[1,1]
        r2 = r2
        
        
        test<- nortest::lillie.test(resid(winning_model, type = "normalized"))
        lilli_test_vector<-test[2][[1]]
        #now test homogenity of variance
        cortests<-cor.test((as.vector(resid(winning_model, type = "normalized"))),
                           as.vector(winning_model$fitted),method = "pearson")
        corfit_vector1<-cortests$p.value
        print("#assign test")
        tables<-data.frame(ID = x,
                           res = pgls_data$res[x],
                           resample_no = d, 
                           p_QWD = p_QWD,
                           p_sdQWD = p_sdQWD,
                           r2 = r2,
                           slope_QWD = slope_QWD, 
                           slope_sdQWD =slope_sdQWD, 
                           slope_error_QWD = slope_error_QWD,
                           slope_error_sdQWD =slope_error_sdQWD,
                           normality = lilli_test_vector,
                           homogeneity = corfit_vector1,
                           intercept =  intercept,
                           exotic_perc = exotic_perc
        )
      }
      
      
      
      
    }
    return(tables)
  })
  
  return(do.call(rbind,resample_pgls))
  
})
resampled_30<-do.call(rbind,results_30)
rownames(resampled_30)<-NULL
saveRDS(resampled_30, "pgls_results/FBD_stem_resampled_30_nonnative_concensus_res4.rds")


resampled_30<-readRDS("pgls_results/FBD_stem_resampled_30_nonnative_concensus_res4.rds")
resampled_25<-readRDS("pgls_results/FBD_stem_resampled_25_nonnative_concensus_res4.rds")
resampled_20<-readRDS("pgls_results/FBD_stem_resampled_20_nonnative_concensus_res4.rds")


source("R/sense_data_function.R")
sd_20<-sense_data(resampled_20, 20)
sd_25<-sense_data(resampled_25, 25)
sd_30<-sense_data(resampled_30, 30)

full_sense<-rbind(sd_20, sd_25, sd_30)
saveRDS(full_sense, "pgls_results/FBD_stem_resampled_nonnative_concensus_res4_sensitivity_data.rds")


