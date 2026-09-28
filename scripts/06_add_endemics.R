
# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 22 June 2026


# Intro ------------------------------------------------------------------------
# This script is used to add the data from:  

# Gallagher, R. V., Allen, S. P., Govaerts, R., Rivers, M. C., Allen, A. P., Keith, D. A., Merow, C., Maitner, B., Butt, N., Auld, T. D., Enquist, B. J., Eiserhardt, W. L.,   Wright, I. J., Mifsud, J. C. O., Espinosa-Ruiz, S., Possingham, H., & Adams, V. M. (2023). Global shortfalls in threat assessments for endemic flora by country. Plants, People, Planet, 5(6), 885–898. https://doi.org/10.1002/ppp3.10369

# It is used to to create estimates of the potential gap that regional Red List assessments for endemic species could fill in the global IUCN Red List.

# Libraries --------------------------------------------------------------------

library(data.table)
library(tidyverse)
library(sf)
library(ggrepel)
library(ggpubr)
library(patchwork)

# Import data ------------------------------------------------------------------

dist_native <- dist_native <- fread("data/dist_native.txt") 

endemics <- dist_native %>% 
  group_by(plant_name_id) %>% 
  filter(n() <2)

plants_full_raw <- fread("data/wcvp_accepted_merged.txt")

growth <- fread("data/growth_forms_nnet.txt") %>% 
  dplyr::select(growth_form, plant_name_id)

redlist_raw <- fread("data/redlist_data_04_2026.csv", sep = ",")

redlist_threat <- redlist_raw %>% 
  filter(redlistCategory %in% c("Endangered", "Vulnerable", "Critically Endangered")) 

tdwg_3 <- st_read(dsn ="data/wgsrpd-master/level3") %>% 
  filter(!LEVEL3_COD == "BOU")

dist_native_end <- dist_native %>% 
  filter(plant_name_id %in% endemics$plant_name_id) %>% 
  filter(area_code_l3 %in% c("CPP"))

xx <- redlist_raw %>% 
  filter(plant_name_id %in% dist_native_end$plant_name_id) 

endemics_toadd <- fread("output/endemics_assessed_05_2026.txt")


threat_search_to_add <- endemics_toadd %>% 
  filter(taxon_name %in% plants_full_raw$taxon_name) %>% 
  filter(!taxon_name %in% redlist_raw$taxon_name) %>% 
  filter(infra_rank == 0) %>% 
  left_join(plants_full_raw, by = "taxon_name") 


#fwrite(threat_search_to_add, "output/redlist_srli/redlist/threat_search_to_add.txt")

threat_search_full <- endemics_toadd %>% 
  filter(taxon_name %in% plants_full_raw$taxon_name) %>% 
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
  group_by(area, area_code_l3) %>% 
  summarise(n_wcvp = n_distinct(plant_name_id))

dist_native_all <- dist_wcvp %>% 
  left_join(dist_threat, by = "area") %>%
  left_join(dist_threat_add, by = "area") %>% 
  
  left_join(dist_red_tot, by = "area") %>% 

  mutate(
         n_comb = (n_red + n_threat_search_add), 
         prop_now  = n_red / n_wcvp, 
         prop_after = n_comb / n_wcvp, 
         prop_diff = prop_after - prop_now, 
         diff = n_comb - n_red)

write.csv(dist_native_all, "redlist_wcvp_threatsearch_comp_table.csv")


dist_native_all_map <- dist_native_all %>% 
  as.data.frame() %>% 
  left_join(tdwg_3, by = c("area_code_l3" = "LEVEL3_COD")) %>% 
  st_as_sf()

  

map_now <- ggplot(dist_native_all_map) +
  geom_sf(aes(fill = (n_red))) +
  paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging", 
                                    breaks = seq(0, 26000, length.out =6), 
                                    limits = c(0, 26000)) +
  labs(fill = "Species", caption = "Global IUCN Red List species only") +
  theme(
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
  ) +
  coord_sf(crs  = sf::st_crs("+proj=eqearth")) 

 map_increase <- ggplot(dist_native_all_map) +
  geom_sf(aes(fill = (n_comb))) +
   paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging", 
                                     breaks = seq(0, 26000, length.out =6), 
                                     limits = c(0, 26000)) +
   labs(fill = "Species", caption ="Global IUCN Red List + Threatsearch endemics") +

   theme(
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
 
 
 # map_diff <- ggplot(dist_native_all_map) +
 #   geom_sf(aes(fill = (prop_diff))) +
 #   #paletteer::scale_fill_paletteer_c("ggthemes::Orange", limits = c(1,25000)) +
 #   #paletteer::scale_fill_paletteer_c("ggthemes::Orange", limits = c(0,1)) +
 #   labs(fill = "Species", caption ="Difference in number of species") +
 #   #paletteer::scale_fill_paletteer_c("pals::ocean.thermal") +
 #   theme_pubclean() +
 #   coord_sf(crs  = sf::st_crs("+proj=eqearth")) +
 #   theme(aspect.ratio = 0.5)
 
 map_tot <- ggplot(dist_native_all_map) +
   paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging", 
                                     breaks = seq(0, 26000, length.out =5), 
                                     limits = c(0, 26000)) +

   geom_sf(aes(fill = (n_wcvp))) +

   labs(fill = "Species", caption ="Total number of species") +

   
   theme(
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

 ggsave("maps_diff_redlist_threatsearch2_lab.png",maps, dpi = 600, width = 14, 
        height = 5)
 


dist_comb <- dist_native %>% 
  filter(plant_name_id %in% threat_search_full$plant_name_id | 
           plant_name_id %in% redlist_raw$plant_name_id) %>% 
  group_by(area_code_l3) %>% 
  summarise(n_assessed = n_distinct(plant_name_id))


dist_miss <- dist_native %>% 
  #filter(plant_name_id %in% threat_search_full$plant_name_id | plant_name_id %in% redlist_raw$plant_name_id) %>% 
  group_by(area_code_l3) %>% 
  summarise(n_wcvp = n_distinct(plant_name_id)) %>% 
  left_join(dist_comb, by = "area_code_l3")  %>% 
  mutate(diff = n_assessed - n_wcvp, 
         prop_miss =  1- n_assessed / n_wcvp) %>% 
  left_join(tdwg_3, by = c("area_code_l3" = "LEVEL3_COD")) %>% 
  st_as_sf() 



prop_missing_total <- ggplot(dist_miss) +
  geom_sf(aes(fill = prop_miss)) +
  paletteer::scale_fill_paletteer_c("grDevices::Plasma")+
  #paletteer::scale_fill_paletteer_c("ggthemes::Orange") +
  #scale_fill_viridis_c() 
  theme_pubclean() +
  coord_sf(crs  = sf::st_crs("+proj=eqearth")) +
  theme(aspect.ratio = 0.5)


ggsave("prop_missing_comb_comp_redlist_threatsearch_pubready.png",
       prop_missing_total, dpi = 600)






threatsearch_names<- plants_full_raw %>% 
  filter(taxon_name %in% threat_search_to_add$taxon_name) 

dist_comb_names  <- dist_native %>% 
  filter(plant_name_id %in% threat_search_full$plant_name_id | plant_name_id %in% redlist_raw$plant_name_id) 


all_growth1 <- growth %>% 
 # filter(plant_name_id %in% redlist_raw$plant_name_id) %>% 
  group_by(growth_form) %>% 
  summarise(Freq = n_distinct(plant_name_id), 
            dataset = "wcvp")

all_growth2 <- growth %>% 
  filter(plant_name_id %in% redlist_raw$plant_name_id) %>% 
  group_by(growth_form) %>% 
  summarise(Freq = n_distinct(plant_name_id), 
            dataset = "redlist")

all_growth3 <- growth %>% 
  filter(plant_name_id %in% redlist_threat$plant_name_id) %>% 
  group_by(growth_form) %>% 
  summarise(Freq = n_distinct(plant_name_id), 
            dataset = "redlist threatened")

all_growth4 <- growth %>% 
  filter(plant_name_id %in% threat_search_to_add$plant_name_id) %>% 
  group_by(growth_form) %>% 
  summarise(Freq = n_distinct(plant_name_id), 
            dataset = "endemics")

all_growth5 <- growth %>% 
  filter(plant_name_id %in% dist_comb_names$plant_name_id) %>% 
  group_by(growth_form) %>% 
  summarise(Freq = n_distinct(plant_name_id), 
            dataset = "endemics + redlist")


all_growth <- rbind(all_growth1, all_growth2, all_growth3, all_growth4, all_growth5)

comb_table <- all_growth %>% 
  group_by(dataset) %>% 
  mutate(prop = Freq / sum(Freq))


wide_table <- pivot_wider(data =  comb_table %>%  dplyr::select(-prop),
                  values_from = c(Freq), names_from = c(growth_form)
                  
                  
                  )

write.csv(wide_table, "output/summary_table_growths_dataset.csv")

