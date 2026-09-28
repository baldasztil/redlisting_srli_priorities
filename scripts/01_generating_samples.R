
# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 22 June 2026


# Intro ------------------------------------------------------------------------
# This script generates random samples of species from the WCVP for 
# the SRLI and the IUCN Red List with the same number of species as the datasets

# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)


# Import data ------------------------------------------------------------------

plants_full_raw <- fread("data/wcvp_accepted_merged.txt")


redlist_raw <- fread("data/redlist_data_04_2026.csv", sep = ",")

srli_raw <-  fread("data/srli_data_09_2026.csv", sep = ",") 
table(srli_raw$redlistCategory)


plants_full <- plants_full_raw 



# manipulate data --------------------------------------------------------------

plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, taxon_name)


redlist_names_raw <- plants_full %>% 
  filter(plant_name_id %in% redlist_raw$plant_name_id) 

srli_names_raw <- plants_full %>% 
  filter(plant_name_id %in% srli_raw$plant_name_id) 


# samples ------------------------------------------------------------------
set.seed(123)
species_samples_red_list <- replicate(10000, sample(x = plants_full$plant_name_id, size = nrow(redlist_names_raw), replace = F), simplify = F)
saveRDS(species_samples_red_list, "output/random_samples/redlist_10000_samples.rds")

species_samples_srli_list <- replicate(10000, sample(x = plants_full$plant_name_id, size = nrow(srli_names_raw), replace = F), simplify = F)
saveRDS(species_samples_srli_list, "output/random_samples/srli_10000_samples.rds")