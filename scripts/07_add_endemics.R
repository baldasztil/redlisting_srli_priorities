library(data.table)
library(tidyverse)
library(sf)
library(ggrepel)

dist_native <- dist_native <- fread("data/dist_native.txt") 

endemics <- dist_native %>% 
  group_by(plant_name_id) %>% 
  filter(n() <2)

plants_full_raw <- fread("data/wcvp_accepted_merged.txt")


redlist_raw <- fread("data/red/cleaned_10_2025/redlist_data_04_2026.csv", sep = ",")
table(redlist_raw$redlistCategory)

srli_raw <-  fread("data/red/cleaned_10_2025/srli_data_04_2026.csv", sep = ",") %>% 
  filter(!is.na(redlistCategory))

length(unique(c(redlist_raw$plant_name_id, srli_raw$plant_name_id)))

srli_threat <-srli_raw %>% 
  filter(redlistCategory %in% c("Endangered", "Vulnerable", "Critically Endangered")) #


redlist_threat <- redlist_raw %>% 
  filter(redlistCategory %in% c("Endangered", "Vulnerable", "Critically Endangered")) # "Endangered", "Vulnerable",


tdwg_3 <- st_read(dsn ="data/wgsrpd-master/level3") %>% 
  filter(!LEVEL3_COD == "BOU")

dist_native_end <- dist_native %>% 
  filter(plant_name_id %in% endemics$plant_name_id) %>% 
  filter(area_code_l3 %in% c("CPP"))


xx <- redlist_raw %>% 
  filter(plant_name_id %in% dist_native_end$plant_name_id) 

xx3 <- fread("output/redlist_srli/redlist/endemics_assessed_05_2026.txt")


threat_search_to_add <- xx3 %>% 
  filter(taxon_name %in% plants_full_raw$taxon_name) %>% 
  filter(!taxon_name %in% redlist_raw$taxon_name) %>% 
  filter(infra_rank == 0) %>% 
  left_join(plants_full_raw, by = "taxon_name") 


fwrite(threat_search_to_add, "output/redlist_srli/redlist/threat_search_to_add.txt")

threat_search_full <- xx3 %>% 
  filter(taxon_name %in% plants_full_raw$taxon_name) %>% 
  #filter(!taxon_name %in% redlist_raw$taxon_name) %>% 
  filter(infra_rank == 0) %>% 
  left_join(plants_full_raw, by = "taxon_name")


dist_threat <- dist_native %>% 
  filter(plant_name_id %in% threat_search_full$plant_name_id) %>% 
  group_by(area) %>% 
  summarise(n_threat_search = n_distinct(plant_name_id))

dist_threat_add <- dist_native %>% 
  filter(plant_name_id %in% threat_search_full$plant_name_id) %>% 
  filter(!plant_name_id %in% redlist_raw$plant_name_id) %>% 
  group_by(area) %>% 
  summarise(n_threat_search_add = n_distinct(plant_name_id))  

dist_red_tot <- dist_native %>% 
    filter(plant_name_id %in% redlist_raw$plant_name_id) %>% 
    group_by(area) %>% 
  summarise(n_red = n_distinct(plant_name_id))

dist_wcvp <- dist_native %>% 
  #filter(plant_name_id %in% threat_search_full$plant_name_id | plant_name_id %in% redlist_raw$plant_name_id) %>% 
  group_by(area, area_code_l3) %>% 
  summarise(n_wcvp = n_distinct(plant_name_id))

dist_native_all <- dist_native %>% 
  filter(plant_name_id %in% endemics$plant_name_id) %>% 
  filter(plant_name_id %in% redlist_raw$plant_name_id) %>% 
  group_by(area) %>% 
  summarise(n_red_end = n_distinct(plant_name_id)) %>% 
  right_join(dist_threat, by = "area") %>%
  right_join(dist_threat_add, by = "area") %>% 
  
  right_join(dist_red_tot, by = "area") %>% 
  right_join(dist_wcvp, by = "area") %>% 
  
  mutate(n_red_end = replace_na(n_red_end, 0), 
         n_comb = (n_red + n_threat_search_add), 
         prop_now  = n_red / n_wcvp, 
         prop_after = n_comb / n_wcvp, 
         prop_diff = prop_after - prop_now, 
         diff = n_comb - n_red)


plot <- ggplot(dist_native_all, aes(x = sqrt(n_threat_search), y = sqrt(n_red_end))) +
 geom_abline(intercept = 0, slope = 1, col = "red", lty = 2) +
  geom_point(aes(col = prop_diff)) +
  geom_text_repel(aes(label = area)) +
  scale_x_continuous(limits = c(0,125), "Number of endemics assessed on global IUCN Red List") +
  scale_y_continuous(limits = c(0,125), "Number of endemics assessed on threat search") +
  theme_bw() +
  labs(color = "Porportional difference") +
  paletteer::scale_color_paletteer_c("ggthemes::Orange", limits = c(0,1)) 
plot

ggsave("endemics_comp_redlist_threatsearch.png",plot, dpi = 600)



plot2 <- ggplot(dist_native_all %>% 
         slice_max(diff, n = 25) %>%  
         arrange(prop_now)%>% 
         mutate(area=factor(area, area))
) +
  
  geom_segment(aes(x= area, xend=area, y=prop_now, yend=prop_after), color="grey") +
  geom_point( aes(x=area, y=prop_now), color=rgb(0.2,0.7,0.1,0.5), size=3 ) +
  geom_point( aes(x=area, y=prop_after), color=rgb(0.7,0.2,0.1,0.5), size=3 ) +
  #geom_abline(intercept = 0, slope = 1, col = "red", lty = 2) +
  #geom_point(aes(col = prop_diff)) +
  #geom_text_repel(aes(label = area)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1, 0.25)) +
  #scale_y_continuous(limits = c(0,1), "Number of endemics assessed on threat search") +
  #labs(col = "Speciesal difference") +
  theme_classic() +
  coord_flip()
    
  

ggsave("lolipop_prop_increase_comp_redlist_threatsearch.png",plot2, dpi = 600)


library(ggpubr)
library(patchwork)
tdwg_3 <- st_read(dsn ="data/wgsrpd-master/level3") %>% 
  filter(!LEVEL3_COD == "BOU")
plot(st_geometry(tdwg_3))

dist_native_all_map <- dist_native_all %>% 
  as.data.frame() %>% 
  left_join(tdwg_3, by = c("area_code_l3" = "LEVEL3_COD")) %>% 
  st_as_sf()


  

map_now <- ggplot(dist_native_all_map) +
  geom_sf(aes(fill = (n_red))) +
  paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging", 
                                    breaks = seq(0, 25000, length.out =6), 
                                    limits = c(0, 25000)) +
  # paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging",
  #                                   breaks = seq(0, 8000, length.out =5 ), 
  #                                   limits = c(0, 8000)) +
 # paletteer::scale_fill_paletteer_c("ggthemes::Orange", limits = c(0,1)) +
  labs(fill = "Species", caption = "Global IUCN Red List species only") +
  
#  paletteer::scale_fill_paletteer_c("grDevices::Sunset") +
 # paletteer::scale_fill_paletteer_c("ggthemes::Orange", limits = c(0,25000)) +
  
  
  theme(
    # panel.spacing = unit(1, "cm"),  # or "cm", "mm", etc.
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(color = "white", face = "bold", size = 6),
    panel.grid.major = element_line(color = "grey85", size = 0.2),
    legend.position = "right",
    legend.box = "vertical",
    axis.text = element_blank(),
    axis.ticks = element_blank(), 
    aspect.ratio = 1/1.5
    # strip.switch.pad.grid = unit(0.4, "cm") 
  ) +
  coord_sf(crs  = sf::st_crs("+proj=eqearth")) 

 map_increase <- ggplot(dist_native_all_map) +
  geom_sf(aes(fill = (n_comb))) +
  # paletteer::scale_fill_paletteer_c("ggthemes::Orange", limits = c(0,25000)) +
   paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging", 
                                     breaks = seq(0, 25000, length.out =6), 
                                     limits = c(0, 25000)) +
   # paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging", 
   #                                   breaks = seq(0, 16000, length.out =5), 
   #                                   limits = c(0, 16000)) +
   
   #paletteer::scale_fill_paletteer_c("ggthemes::Orange", limits = c(0,1)) +
   labs(fill = "Species", caption ="Global IUCN Red List + Threatsearch endemics") +
   
  # paletteer::scale_fill_paletteer_c("grDevices::Sunset") +
  # paletteer::scale_fill_paletteer_c("pals::ocean.thermal") +
   
   theme(
     # panel.spacing = unit(1, "cm"),  # or "cm", "mm", etc.
     panel.background = element_rect(fill = "white", color = NA),
     plot.background = element_rect(fill = "white", color = NA),
     strip.background = element_rect(fill = "black"),
     strip.text = element_text(color = "white", face = "bold", size = 6),
     panel.grid.major = element_line(color = "grey85", size = 0.2),
     legend.position = "right",
     legend.box = "vertical",
     axis.text = element_blank(),
     axis.ticks = element_blank(), 
     aspect.ratio = 1/1.5
     # strip.switch.pad.grid = unit(0.4, "cm") 
   ) +
  coord_sf(crs  = sf::st_crs("+proj=eqearth"))
 # 
 # map_diff <- ggplot(dist_native_all_map) +
 #   geom_sf(aes(fill = (prop_diff))) +
 #   #paletteer::scale_fill_paletteer_c("ggthemes::Orange", limits = c(1,25000)) +
 #   #paletteer::scale_fill_paletteer_c("ggthemes::Orange", limits = c(0,1)) +
 #   labs(fill = "Species", caption ="Difference in number of species") +
 #   #paletteer::scale_fill_paletteer_c("pals::ocean.thermal") +
 #   theme_pubclean() +
 #   coord_sf(crs  = sf::st_crs("+proj=eqearth")) +
 #   theme(aspect.ratio = 0.5)
 # 
 map_tot <- ggplot(dist_native_all_map) +
   paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging", 
                                     breaks = seq(0, 26000, length.out =5), 
                                     limits = c(0, 26000)) +

   geom_sf(aes(fill = (n_wcvp))) +
   # paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging", direction = -1, 
   #                                   limits = c(1, 25000)) +  
   #paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging") 
  # paletteer::scale_fill_paletteer_c("ggthemes::Orange", limits = c(0,25000)) +
   labs(fill = "Species", caption ="Total number of species") +
   #paletteer::scale_fill_paletteer_c("grDevices::Sunset") +
 #  paletteer::scale_fill_paletteer_c("pals::ocean.thermal") +
   
   
   theme(
     # panel.spacing = unit(1, "cm"),  # or "cm", "mm", etc.
     panel.background = element_rect(fill = "white", color = NA),
     plot.background = element_rect(fill = "white", color = NA),
     strip.background = element_rect(fill = "black"),
     strip.text = element_text(color = "white", face = "bold", size = 6),
     panel.grid.major = element_line(color = "grey85", size = 0.2),
     legend.position = "right",
     legend.box = "vertical",
     axis.text = element_blank(),
     axis.ticks = element_blank(), 
     aspect.ratio = 1/1.5
     # strip.switch.pad.grid = unit(0.4, "cm") 
   ) +
   coord_sf(crs  = sf::st_crs("+proj=eqearth")) 

 
 maps <- map_now + map_increase + map_tot  + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
   theme(plot.tag = element_text(face = 'bold'))

 ggsave("maps_diff_redlist_threatsearch2_lab.png",maps, dpi = 600, width = 14, height = 5)
 


dist_comb <- dist_native %>% 
  filter(plant_name_id %in% threat_search_full$plant_name_id | plant_name_id %in% redlist_raw$plant_name_id) %>% 
  group_by(area_code_l3) %>% 
  summarise(n_assessed = n_distinct(plant_name_id))




dist_wcvp <- dist_native %>% 
  #filter(plant_name_id %in% threat_search_full$plant_name_id | plant_name_id %in% redlist_raw$plant_name_id) %>% 
  group_by(area_code_l3) %>% 
  summarise(n_wcvp = n_distinct(plant_name_id)) %>% 
  left_join(dist_comb, by = "area_code_l3")  %>% 
  mutate(diff = n_assessed - n_wcvp, 
         prop_miss =  1- n_assessed / n_wcvp) %>% 
  left_join(tdwg_3, by = c("area_code_l3" = "LEVEL3_COD")) %>% 
  st_as_sf() 



prop_missing_total <- ggplot(dist_wcvp) +
  geom_sf(aes(fill = prop_miss)) +
  #paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging") 
  paletteer::scale_fill_paletteer_c("ggthemes::Orange") +
  #scale_fill_viridis_c() 
  theme_pubclean() +
  coord_sf(crs  = sf::st_crs("+proj=eqearth")) +
  theme(aspect.ratio = 0.5)


ggsave("prop_missing_comb_comp_redlist_threatsearch.png",prop_missing_total, dpi = 600)






threatsearch_names<- plants_full %>% 
  filter(taxon_name %in% threat_search_to_add$taxon_name) 


all_growth <- as.data.frame(growth.form.rand(species_sample =  threatsearch_names$plant_name_id, method = "overall"))

all_growth$pro_freq  <- all_growth$Freq / sum(all_growth$Freq)

names(all_growth)[1] <- "growth_form"


dist_comb_names  <- dist_native %>% 
  filter(plant_name_id %in% threat_search_full$plant_name_id | plant_name_id %in% redlist_raw$plant_name_id) 

all_growth <- as.data.frame(growth.form.rand(species_sample =  dist_comb_names$plant_name_id, method = "overall"))

all_growth$pro_freq  <- all_growth$Freq / sum(all_growth$Freq)

names(all_growth)[1] <- "growth_form"


all_growth <- as.data.frame(growth.form.rand(species_sample =  dist_native$plant_name_id, method = "overall"))

all_growth$pro_freq  <- all_growth$Freq / sum(all_growth$Freq)

names(all_growth)[1] <- "growth_form"



all_growth <- as.data.frame(growth.form.rand(species_sample =  redlist_raw$plant_name_id, method = "overall"))

all_growth$pro_freq  <- all_growth$Freq / sum(all_growth$Freq)

names(all_growth)[1] <- "growth_form"

