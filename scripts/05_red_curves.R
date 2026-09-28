
# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 22 June 2026


# Intro ------------------------------------------------------------------------
# This script is used to create accumulation curves of red list assessments
# through time

# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(sf)
library(terra)
options(dplyr.summarise.inform = FALSE) 


# Import data ------------------------------------------------------------------

dist_native <- dist_native <- fread("data/dist_native.txt") 
plants_full_raw <- fread("data/wcvp_accepted_merged.txt")


redlist_raw <- fread("data/redlist_data_04_2026.csv", sep = ",")
table(redlist_raw$redlistCategory)




tdwg_3 <- st_read(dsn ="data/wgsrpd-master/level3") %>% 
  filter(!LEVEL3_COD == "BOU")




growth <- fread("data/growth_forms_nnet.txt") %>% 
  dplyr::select(growth_form, plant_name_id)

plants_full <- plants_full_raw %>% 
  left_join(growth, by = "plant_name_id")



lookup <- c(LEVEL2_NAM = "region.x", LEVEL1_NAM = "continent.x")

dist_native_calc <- dist_native %>% 
  left_join(plants_full, by = "plant_name_id") %>% 
  rename(all_of(lookup))

prop_org <- dist_native_calc %>%  
  group_by(LEVEL2_NAM, LEVEL1_NAM) %>% 
  summarise(prop_wcvp = n_distinct(plant_name_id) / nrow(plants_full))



# manipulate data --------------------------------------------------------------

plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, taxon_name, growth_form)


plantlist_dist_phylo_growth <- dist_native %>% 
  left_join(plantlist_names, by = "plant_name_id") %>% 
  dplyr::select(plant_name_id, taxon_name, LEVEL3_COD = area_code_l3, growth_form, 
                family,LEVEL1_NAM = continent, 
                LEVEL2_NAM = region) %>% 
  mutate(label = gsub(" ", "_", taxon_name))



# manipulate data --------------------------------------------------------------
options(dplyr.summarise.inform = FALSE) 



all_growth <-plants_full %>% 
  group_by(growth_form) %>% 
  summarise(Freq = n())


# curves_g <- redlist_raw %>% 
#   left_join(growth, by = "plant_name_id") %>% 
#   group_by(yearPublished, growth_form) %>% 
#   summarise(n = n_distinct(plant_name_id), .groups = "drop") %>% 
#   drop_na(.) %>% 
#   complete(
#     yearPublished = full_seq(yearPublished, 1),      # fills all yearPublisheds in sequence
#     growth_form,                   # ensures each growth form appears every yearPublished
#     fill = list(n = 0)             # replace missing with 0
#   ) %>% 
#   arrange(growth_form, yearPublished) %>% 
#   group_by(growth_form) %>% 
#   mutate(sum = cumsum(n)) %>% 
#   left_join(all_growth, by = c("growth_form")) %>% 
#   mutate(prop_sum = cumsum(n) / Freq)  


div <- dist_native_calc %>% 
  group_by(LEVEL1_NAM) %>% 
  summarise(div = n_distinct(plant_name_id), .groups = "drop")

curves_div <- dist_native_calc %>% 
  left_join(redlist_raw, by = "plant_name_id") %>% 
  # left_join(growth, by = "plant_name_id") %>% 
  group_by(yearPublished, LEVEL1_NAM) %>% 
  summarise(n = n_distinct(plant_name_id), .groups = "drop") %>% 
  drop_na(.) %>% 
  complete(
    yearPublished = full_seq(yearPublished, 1),      # fills all yearPublisheds in sequence
    LEVEL1_NAM,                   # ensures each growth form appears every yearPublished
    fill = list(n = 0)             # replace missing with 0
  ) %>% 
  arrange(LEVEL1_NAM, yearPublished) %>% 
  group_by(LEVEL1_NAM) %>% 
  mutate(sum = cumsum(n))  %>% 
  left_join(div, by = c("LEVEL1_NAM")) %>% 
  mutate(prop_sum = cumsum(n) / div)  

fwrite(curves_div, "red_continent_acc_curves_09_2026.txt")
#fwrite(curves_g, "red_growth_acc_curves_04_2026.txt")

 