library(tidyverse)
library(patchwork)
library(h3jsr)
library(sf)
#load data

differences4<-readRDS("difference_data/differences_res4.rds")
differences3<-readRDS("difference_data/differences_res3.rds")
differences5<-readRDS("difference_data/differences_res5.rds")


states_poly<-st_read("geo_data/States_shapefile-shp/States_shapefile.shp")
states_poly<-states_poly %>%
  st_as_sf() %>%
  st_transform(crs = "EPSG:5070")



differences4=lapply(differences4, function(x){
  x<-x%>% mutate(resolution = "~1,770km2")  
  poly<-cell_to_polygon(x$res)
  poly<-cbind(x,poly) %>%
    st_as_sf() %>%
    st_transform(crs = "EPSG:5070")
    
  return(poly)
})
differences5=lapply(differences5, function(x){
  x<-x%>% mutate(resolution = "~250km2")  
  poly<-cell_to_polygon(x$res)
  poly<-cbind(x,poly) %>%
    st_as_sf() %>%
    st_transform(crs = "EPSG:5070")
  
  return(poly)
})
differences3=lapply(differences3, function(x){
  x<-x%>% mutate(resolution = "~12,330km2")  
  poly<-cell_to_polygon(x$res)
  poly<-cbind(x,poly) %>%
    st_as_sf() %>%
    st_transform(crs = "EPSG:5070")
  
  return(poly)
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

diff_total<-rbind(differences4,differences5,differences3)
diff_total<-diff_total %>%
  mutate(PGLS_sample_size = as.factor(PGLS_sample_size)) %>%
  mutate(name = fct_relevel(resolution, 
                            "~250km2", "~1,770km2", "~12,330km2")) 



L = st_cast(states_poly,"MULTILINESTRING")





levels(diff_total$name)<-c(expression("250"~km^2),
                           expression("1,770"~km^2),
                           expression("12,330"~km^2))
#total

total<-ggplot() +
  geom_sf(data = states_poly, fill = "black", color = "white") +
  geom_sf(data =diff_total, 
          aes(fill = slope_diff, color = slope_diff)) +
  geom_sf(data = L,color = "white") +
  scale_x_continuous(breaks = c(-88, -85, -82)) +
  facet_grid(PGLS_sample_size ~ name,
             labeller = label_parsed) +
  theme_minimal() +
  labs(fill = "Slope\ndifference",
       color = "Slope\ndifference") +
  scale_fill_gradient2(
    low = "deepskyblue", mid = "azure3", high = "orange",
    midpoint = 0,
    name = "Slope\ndifference"
  ) + 
  scale_color_gradient2(
    low = "deepskyblue", mid = "azure3", high = "orange",
    midpoint = 0,
    name = "Slope\ndifference"
  ) + 
  coord_sf(xlim = c(605755.4,
                    1632548.473),
           ylim = c(
             1406449.988, 250506.46)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey",  linetype = 3),
        legend.position = "top",
        legend.background = element_blank(),
        legend.title = element_text(face = "bold", size = 14),
        legend.direction = "horizontal",
        plot.background = element_blank()) 

total



total_non<-ggplot() +
  geom_sf(data = states_poly, fill = "black", color = "white") +
  geom_sf(data =diff_total, 
          aes(fill = exotic_perc_avg *100, 
              color = exotic_perc_avg *100)) +
  geom_sf(data = L,color = "white") +
  scale_x_continuous(breaks = c(-88, -85, -82)) +
  facet_grid(PGLS_sample_size ~ name,
             labeller = label_parsed) +
  theme_minimal() + 
  labs(fill = "Nonnative\nant %", color = "Nonnative\nant %") +
 scale_fill_viridis_c(option = "H") + 
  scale_color_viridis_c(option = "H") + 
  coord_sf(xlim = c(605755.4,
                    1632548.473),
           ylim = c(
             1406449.988, 250506.46)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey",  linetype = 3),
        legend.position = "top",
        legend.title = element_text(face = "bold", size = 14),
        legend.direction = "horizontal",
        plot.background = element_blank()) 

total_non

total + total_non + plot_layout(nrow = 1) &
  theme(legend.box.spacing = unit(0, "pt"),
        plot.background = element_blank())


##sep
small<-ggplot() +
  geom_sf(data = states_poly, fill = "black", color = "white") +
  geom_sf(data =differences5, 
          aes(fill = slope_diff, color = slope_diff)) +
  geom_sf(data = L,color = "white") +
  facet_wrap(vars(PGLS_sample_size), nrow = 3) +
  theme_bw() +
  scale_fill_gradient2(
    low = "deepskyblue", mid = "azure3", high = "orange",
    midpoint = 0,
    name = "Slope\ndifference"
  ) + 
  scale_color_gradient2(
    low = "deepskyblue", mid = "azure3", high = "orange",
    midpoint = 0,
    name = "Slope\ndifference"
  ) + 
  coord_sf(xlim = c(605755.4,
                    1632548.473),
           ylim = c(
             1406449.988, 250506.46)) +
  labs(title = "Spatial resolution:\n~ 252km2") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey",  linetype = 3),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.title = element_text(face = "bold", size = 14),
        legend.direction = "vertical") 
  

med<-ggplot() +
  geom_sf(data = states_poly, fill = "black", color = "white") +
  geom_sf(data =differences4, 
          aes(fill = slope_diff, color = slope_diff)) +
  geom_sf(data = L,color = "white") +
  facet_wrap(vars(PGLS_sample_size),
             nrow = 3) +
  theme_bw() +
  scale_fill_gradient2(
    low = "deepskyblue", mid = "azure3", high = "orange",
    midpoint = 0,
    name = "Slope\ndifference"
  ) + 
  scale_color_gradient2(
    low = "deepskyblue", mid = "azure3", high = "orange",
    midpoint = 0,
    name = "Slope\ndifference"
  ) + 
  coord_sf(xlim = c(605755.4,
                    1632548.473),
           ylim = c(
             1406449.988, 250506.46)) +
  labs(title = "Spatial resolution:\n~ 1,770km2") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey",  linetype = 3),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.title = element_text(face = "bold", size = 14),
        legend.direction = "vertical") 


large<-ggplot() +
  geom_sf(data = states_poly, fill = "black", color = "white") +
  geom_sf(data =differences3, 
          aes(fill = slope_diff, color = slope_diff)) +
  geom_sf(data = L,color = "white") +
  facet_wrap(vars(PGLS_sample_size), nrow = 3) +
  theme_bw() +
  scale_fill_gradient2(
    low = "deepskyblue", mid = "azure3", high = "orange",
    midpoint = 0,
    name = "Slope\ndifference"
  ) + 
  scale_color_gradient2(
    low = "deepskyblue", mid = "azure3", high = "orange",
    midpoint = 0,
    name = "Slope\ndifference"
  ) + 
  coord_sf(xlim = c(605755.4,
                    1632548.473),
           ylim = c(
             1406449.988, 250506.46)) +
  labs(title = "Spatial resolution:\n~ 12,330km2") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey",  linetype = 3),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.title = element_text(face = "bold", size = 14),
        legend.direction = "vertical") 
small +med + large + plot_layout(guides = "collect")

