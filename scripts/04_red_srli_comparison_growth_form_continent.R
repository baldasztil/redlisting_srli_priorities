# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 22 June 2026


# Intro ------------------------------------------------------------------------
# This script is used to compare the growth form distribution of a random samples 
# of species from the WCVP to the IUCN Red List and the SRLI 

# Intro ------------------------------------------------------------------------


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(sf)
library(spmodel)
library(ape)
library(phyloregion)
library(GWmodel)
library(feather)
library(vegan)
library(gmodels)
library(hillR)
library(ggpmisc)
library(RColorBrewer)
library(colorBlindness)
library(tmap)
library(moments)
library(paletteer)
library(factoextra)
library(FSA)
library(FactoMineR)
library(lemon)
library(ggdensity)
library(ggpointdensity)
library(ggblend)
library(geomtextpath)
library(ggrepel)
library(ggdist)
library(ggridges)
library(rstatix)
library(plotrix)
library(ggpubr)
library(egg)
library(patchwork)

growth.form.rand <- function(species_sample, method) {
  #species_sample <- sample_n(plantlist_names, spec_n)
  
  dist <- plantlist_dist_phylo_growth %>% 
    filter(plant_name_id %in% species_sample) 
  
  table_data_sample_random_pre <- dist %>% 
    group_by(plant_name_id, growth_form, LEVEL1_NAM) %>% 
    summarise()
  
  table_data_sample_random <- table(table_data_sample_random_pre$LEVEL1_NAM,
                                    table_data_sample_random_pre$growth_form)
  
  growth_random <- as.data.frame.matrix(table_data_sample_random) 
  
  
  growth_random <- growth_random %>% 
    mutate(LEVEL1_NAM = rownames(growth_random)) 
  
  dist_sum <- dist %>% 
    group_by(plant_name_id, growth_form) %>% 
    summarise()
  table_overall_sample_random <- table(dist_sum$growth_form)
  
  out <- list(table_overall_sample_random, growth_random)
  
  if(method == "overall") {
    return(table_overall_sample_random)
  }
  
  if(method == "continent") {
    return(growth_random)
  }
  
  if(is.null(method)) {
    return(out)
  }
  
}


# Import data ------------------------------------------------------------------

dist_native <- dist_native <- fread("data/dist_native.txt") 
plants_full_raw <- fread("data/wcvp_accepted_merged.txt")


redlist_raw <- fread("data/redlist_data_04_2026.csv", sep = ",")
table(redlist_raw$redlistCategory)

srli_raw <-  fread("data/srli_data_04_2026.csv", sep = ",") %>% 
  filter(!is.na(redlistCategory))

length(unique(c(redlist_raw$plant_name_id, srli_raw$plant_name_id)))


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

plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, taxon_name, growth_form)


plantlist_dist_phylo_growth <- dist_native %>% 
  left_join(plantlist_names, by = "plant_name_id") %>% 
  dplyr::select(plant_name_id, taxon_name, LEVEL3_COD = area_code_l3, growth_form, family,LEVEL1_NAM = continent, 
                LEVEL2_NAM = region) %>% 
  mutate(label = gsub(" ", "_", taxon_name))


redlist_names_raw <- plants_full %>% 
  filter(plant_name_id %in% redlist_raw$plant_name_id)


srli_names_raw <- plants_full %>% 
  filter(plant_name_id %in% srli_raw$plant_name_id) 


species_samples_red_list <- readRDS("output/random_samples/redlist_10000_samples.rds")
species_samples_srli_list <- readRDS("output/random_samples/srli_10000_samples.rds")



# manipulate data --------------------------------------------------------------
options(dplyr.summarise.inform = FALSE) 

all_growth <- as.data.frame(growth.form.rand(species_sample =  plants_full$plant_name_id, method = "overall"))
names(all_growth)[1] <- "growth_form"


rand_growth <- as.data.frame(growth.form.rand(species_sample = redlist_names_raw$plant_name_id,method = "overall"))
names(rand_growth)[1] <- "growth_form"



redlist_null_all <- lapply(species_samples_red_list, growth.form.rand, 
                           method = "overall") 

red_data_all <- do.call(bind_rows, redlist_null_all) %>% 
  mutate(dataset = "redlist")

rand_growth_f <- do.call(bind_rows, redlist_null_all) %>% 
  pivot_longer(cols = names(.), names_to = "growth_form") %>%
  as.data.frame() %>% 
  left_join(rand_growth, by = c("growth_form")) %>% 
  group_by(growth_form) %>% 
  summarise(mean = mean(value), 
            prop_mean = mean(value / Freq), 
            
            q_lower = quantile(value, probs = 0.025), 
            q_upper = quantile(value, probs = 0.975), 
            prop_q_lower = quantile(value / Freq, probs = 0.025), 
            prop_q_upper = quantile(value / Freq, probs = 0.975)) 


red_growth <- plantlist_dist_phylo_growth %>% 
  filter(plant_name_id %in% redlist_names_raw$plant_name_id) %>% 
  left_join(all_growth, by = c("growth_form")) %>% 
  group_by(growth_form) %>% 
  summarise(n = n_distinct(plant_name_id), 
            n_prop = n / unique(Freq))


dif_to_baseline_red <- red_growth %>% 
  left_join(rand_growth_f, by = c("growth_form")) %>% 
  mutate(
    dif = n - mean, 
    dif_max = n - q_lower,
    dif_min = n - q_upper,
    
    dif2 = n_prop - prop_mean, 
    dif2_max = n_prop - prop_q_lower,
    dif2_min = n_prop - prop_q_upper,
    
    dif_prop = (n_prop / prop_mean) - 1, 
    dif_prop_max = (n_prop / prop_q_lower)-1,
    dif_prop_min = (n_prop / prop_q_upper) -1,
    dataset = "redlist"
    
  ) 




# continental patterns  --------------------------------------------------------
options(dplyr.summarise.inform = FALSE)

rand_growth <- growth.form.rand(species_sample = (redlist_names_raw$plant_name_id),
                                method = "continent")


redlist_null_cont <- lapply(species_samples_red_list, growth.form.rand, 
                            method = "continent") 

red_data_cont <- do.call(bind_rows, redlist_null_cont) %>% 
  mutate(dataset = "redlist")

growths <-unique(growth$growth_form)
rand_growth_f <- do.call(bind_rows, redlist_null_cont) %>% 
  pivot_longer(cols = all_of(growths), names_to = "growth_form") %>%
  as.data.frame() %>% 
  group_by(growth_form, LEVEL1_NAM) %>% 
  summarise(mean = mean(value), 
            sd = sd(value), 
            
            min_sd = mean - sd, 
            max_sd = mean + sd, 
            
            q_lower = quantile(value, probs = 0.025), 
            q_upper = quantile(value, probs = 0.975))

all_growth_cont <- as.data.frame(growth.form.rand(species_sample = plants_full$plant_name_id,method = "continent")) %>% 
  pivot_longer(growths, names_to = "growth_form")

redlist_continent_growths <- plantlist_dist_phylo_growth %>% 
  filter(plant_name_id %in% redlist_names_raw$plant_name_id) %>% 
  group_by(growth_form, LEVEL1_NAM) %>% 
  summarise(n = n_distinct(plant_name_id)) %>% 
  right_join(rand_growth_f, by = c("LEVEL1_NAM", "growth_form")) %>% 
  left_join(all_growth_cont, by = c("growth_form", "LEVEL1_NAM")) %>% 
  mutate(
    
    dif_total = value - n, 
    dif_total_prop = (value - n) / value, 
    
    
    dif = n - mean, 
    dif_max = n - q_lower,
    dif_min = n - q_upper,
    
    dif2 =  (n - mean) / mean, 
    dif2_max = (n - q_lower) /  q_lower,
    dif2_min = (n - q_upper) / q_upper,
    
    dif3 =  (n - mean) / sd, 
    dif3_max = (n - q_lower) /  sd,
    dif3_min = (n - q_upper) / sd,
    
    
    dif_prop = (n / nrow(srli_names_raw)) - (mean / nrow(srli_names_raw)),
    dif_prop_max = (n / nrow(srli_names_raw)) - (q_lower / nrow(srli_names_raw)),
    dif_prop_min = (n / nrow(srli_names_raw)) - (q_upper / nrow(srli_names_raw)), 
    
    dataset = "redlist"
  ) %>% 
  mutate(sig = !(dif_min <= 0 & dif_max >= 0))



# srli patterns  ---------------------------------------------------------------



srli_null_all <- lapply(species_samples_srli_list, growth.form.rand,  
                        method = "overall") 

red_srli_data_all <- do.call(bind_rows, srli_null_all) %>% 
  mutate(dataset = "srli") %>% 
  rbind(red_data_all)


rand_growth_f <- do.call(bind_rows, srli_null_all) %>% 
  pivot_longer(cols = names(.), names_to = "growth_form") %>%
  as.data.frame() %>% 
  left_join(all_growth, by = c("growth_form")) %>% 
  group_by(growth_form) %>% 
  summarise(mean = mean(value), 
            prop_mean = mean(value / Freq), 
            
            q_lower = quantile(value, probs = 0.025), 
            q_upper = quantile(value, probs = 0.975), 
            prop_q_lower = quantile(value / Freq, probs = 0.025), 
            prop_q_upper = quantile(value / Freq, probs = 0.975)) 


srli_growth <- plantlist_dist_phylo_growth %>% 
  filter(plant_name_id %in% srli_names_raw$plant_name_id) %>% 
  left_join(all_growth, by = c("growth_form")) %>% 
  group_by(growth_form) %>% 
  summarise(n = n_distinct(plant_name_id), 
            n_prop = n / unique(Freq))


dif_to_baseline_srli <- srli_growth %>% 
  left_join(rand_growth_f, by = c("growth_form")) %>% 
  mutate(
    dif = n - mean, 
    dif_max = n - q_lower,
    dif_min = n - q_upper,
    
    dif2 = n_prop - prop_mean, 
    dif2_max = n_prop - prop_q_lower,
    dif2_min = n_prop - prop_q_upper,
    
    dif_prop = (n_prop / prop_mean) - 1, 
    dif_prop_max = (n_prop / prop_q_lower)-1,
    dif_prop_min = (n_prop / prop_q_upper) -1,
    dataset = "srli"
    
  ) 



rand_growth <- growth.form.rand(species_sample = (srli_names_raw$plant_name_id),
                                method = "continent")


srli_null_cont <- lapply(species_samples_srli_list, growth.form.rand, 
                         method = "continent") 

red_srli_data_cont <- do.call(bind_rows, srli_null_cont) %>% 
  mutate(dataset = "srli") %>% 
  rbind(red_data_cont)


# srli_data <- do.call(bind_rows, srli_null_cont)
# fwrite(srli_data, "red_srli_null_controwth.txt")

growths <-unique(growth$growth_form)
rand_growth_f <- do.call(bind_rows, srli_null_cont) %>% 
  pivot_longer(cols = all_of(growths), names_to = "growth_form") %>%
  as.data.frame() %>% 
  group_by(growth_form, LEVEL1_NAM) %>% 
  summarise(mean = mean(value), 
            sd = sd(value), 
            
            min_sd = mean - sd, 
            max_sd = mean + sd, 
            
            q_lower = quantile(value, probs = 0.025), 
            q_upper = quantile(value, probs = 0.975))

all_growth_cont <- as.data.frame(growth.form.rand(species_sample = plants_full$plant_name_id,method = "continent")) %>% 
  pivot_longer(growths, names_to = "growth_form")

srli_continent_growths <- plantlist_dist_phylo_growth %>% 
  filter(plant_name_id %in% srli_names_raw$plant_name_id) %>% 
  group_by(growth_form, LEVEL1_NAM) %>% 
  summarise(n = n_distinct(plant_name_id)) %>% 
  right_join(rand_growth_f, by = c("LEVEL1_NAM", "growth_form")) %>% 
  left_join(all_growth_cont, by = c("growth_form", "LEVEL1_NAM")) %>% 
  mutate(
    
    dif_total = value - n, 
    dif_total_prop = (value - n) / value, 
    
    
    dif = n - mean, 
    dif_max = n - q_lower,
    dif_min = n - q_upper,
    
    dif2 =  (n - mean) / mean, 
    dif2_max = (n - q_lower) /  q_lower,
    dif2_min = (n - q_upper) / q_upper,
    
    dif3 =  (n - mean) / sd, 
    dif3_max = (n - q_lower) /  sd,
    dif3_min = (n - q_upper) / sd,
    
    
    dif_prop = (n / nrow(srli_names_raw)) - (mean / nrow(srli_names_raw)),
    
    dif_prop_max = (n / nrow(srli_names_raw)) - (q_lower / nrow(srli_names_raw)),
    dif_prop_min = (n / nrow(srli_names_raw)) - (q_upper / nrow(srli_names_raw)), 
    
    dataset = "srli"
    
    
  ) %>% 
  mutate(sig = !(dif_min <= 0 & dif_max >= 0))



dif_to_baseline <- rbind(dif_to_baseline_red, dif_to_baseline_srli)
continent_growths <- rbind(redlist_continent_growths, srli_continent_growths)


fwrite(red_srli_data_cont, "red_srli_null_growth_04_2026.txt")
fwrite(red_srli_data_cont, "red_srli_null_growth_cont_04_2026.txt")

fwrite(dif_to_baseline, "red_srli_null_growth_stats_04_2026.txt")
fwrite(continent_growths, "red_srli_null_growth_stats_cont_04_2026.txt")


