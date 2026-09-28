
# Authors: Ludwig Baldaszti - lbaldaszti@rbge.org.uk
# Date: 22 June 2026


# Intro ------------------------------------------------------------------------
# This script is used to compare the properties of a random samples of
# species from the WCVP to the IUCN Red List and the SRLI 

# Libraries --------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(sf)
library(terra)
options(dplyr.summarise.inform = FALSE) 


subsampling.plants.null <- function(source, prop_table) {
  
  species_sample <- sample_n(plants_full, nrow(source), replace = F)
  
  
  dist <- plantlist_dist_phylo_growth %>% 
    filter(plant_name_id %in% species_sample$plant_name_id)
  
  sp_num <- dist %>%
    right_join(continent_names, by = c("LEVEL3_COD", "LEVEL1_NAM")) %>% 
    group_by(LEVEL3_COD, LEVEL1_NAM) %>%  
    summarise(n = n_distinct(plant_name_id)) %>% 
    mutate(n = ifelse(is.na(n), 0, n), 
           n_prop = n / nrow(source)) %>%  
  left_join(prop_table, by = c("LEVEL3_COD", "LEVEL1_NAM")) %>% 
  mutate(prop_diff = (unique(prop_org) - n_prop))
  
  
  return(sp_num)
}


subsampling.plants.null.cont <- function(species_sample, source) {
  
  #species_sample <- sample_n(plants_full, nrow(source), replace = F)
  

  dist <- plantlist_dist_phylo_growth %>% 
    filter(plant_name_id %in% species_sample)
  
  
  prop_table <- plantlist_dist_phylo_growth %>% 
    filter(plant_name_id %in% source$plant_name_id) %>% 
    group_by(LEVEL1_NAM) %>%  
    summarise(prop_org = n_distinct(plant_name_id) / nrow(source)) %>% 
    ungroup() %>% 
    mutate(prop_org = ifelse(is.na(prop_org), 0, prop_org), 
           prop_org_std = prop_org / sum(prop_org))
  
  
  sp_num <- dist %>%
    group_by(LEVEL1_NAM) %>%  
    summarise(n = n_distinct(plant_name_id)) %>% 
    mutate(n = ifelse(is.na(n), 0, n), 
           n_prop = n / nrow(source), 
           n_prop_std = n_prop / sum(n_prop)) %>%  
    left_join(prop_table, by = c("LEVEL1_NAM")) %>% 
    mutate(prop_diff = (unique(prop_org) - n_prop), 
           prop_diff_std = (unique(prop_org_std) - n_prop_std))
  
  
  return(sp_num)
}

# Import data ------------------------------------------------------------------

dist_native <- dist_native <- fread("data/dist_native.txt") 
plants_full_raw <- fread("data/wcvp_accepted_merged.txt")


redlist_raw <- fread("data/redlist_data_04_2026.csv", sep = ",")

redlist_threat <- redlist_raw %>% 
  filter(redlistCategory %in% c("Endangered", "Vulnerable", "Critically Endangered")) # "Endangered", "Vulnerable",

redlist_ext <- redlist_raw %>% 
  filter(redlistCategory %in% c("Extinct", "Extinct in the Wild"))

redlist_ntlc <- redlist_raw %>%
  filter(redlistCategory %in% c("Near Threatened", "Least Concern"))

redlist_dd <- redlist_raw %>% 
  filter(redlistCategory %in% c("Data Deficient")) #



srli_raw <-  fread("data/srli_data_09_2026.csv", sep = ",") 


srli_ext <- srli_raw %>% 
  filter(redlistCategory %in% c("Extinct", "Extinct in the Wild"))

srli_threat <-srli_raw %>% 
  filter(redlistCategory %in% c("Endangered", "Vulnerable", "Critically Endangered")) #

srli_dd <-srli_raw %>% 
  filter(redlistCategory %in% c("Data Deficient"))  

tdwg_3 <- st_read(dsn ="data/wgsrpd-master/level3") %>% 
  filter(!LEVEL3_COD == "BOU")


xxx <- srli_raw %>% 
  filter(year > 2015)
  

continent_names <- tdwg_3 %>% 
  dplyr::select(LEVEL1_NAM, LEVEL3_COD, LEVEL2_NAM) %>% 
  st_drop_geometry() %>% 
  dplyr::group_by(LEVEL1_NAM, LEVEL3_COD, LEVEL2_NAM) %>% 
  summarise()
  

growth <- fread("data/growth_forms_nnet.txt") %>% 
  dplyr::select(growth_form, plant_name_id)

plants_full <- plants_full_raw %>% 
  left_join(growth, by = "plant_name_id")



lookup <- c(LEVEL3_COD = "area_code_l3", LEVEL1_NAM = "continent.x")

dist_native_calc <- dist_native %>% 
  left_join(plants_full, by = "plant_name_id") %>% 
  rename(all_of(lookup))

prop_org <- dist_native_calc %>%  
  group_by(LEVEL3_COD, LEVEL1_NAM) %>% 
  summarise(prop_wcvp = n_distinct(plant_name_id) / nrow(plants_full))
  


# manipulate data --------------------------------------------------------------

plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, taxon_name, growth_form)


plantlist_dist_phylo_growth <- dist_native %>% 
  left_join(plantlist_names, by = "plant_name_id") %>% 
  dplyr::select(plant_name_id, taxon_name, LEVEL3_COD = area_code_l3, growth_form, 
                family,LEVEL1_NAM = continent) %>% 
  mutate(label = gsub(" ", "_", taxon_name))



redlist_names_raw <- plants_full %>% 
  filter(plant_name_id %in% redlist_raw$plant_name_id) 

redlist_names_threat <- plants_full %>% 
  filter(plant_name_id %in% redlist_threat$plant_name_id) 


prop_redlist <- plantlist_dist_phylo_growth %>% 
  filter(plant_name_id %in% redlist_names_raw$plant_name_id) %>% 
  group_by(LEVEL3_COD, LEVEL1_NAM) %>%  
  summarise(prop_org = n_distinct(plant_name_id) / nrow(redlist_names_raw)) %>% 
  right_join(continent_names, by = c("LEVEL3_COD", "LEVEL1_NAM"))  %>% 
  mutate(prop_org = ifelse(is.na(prop_org), 0, prop_org))


srli_names_raw <- plants_full %>% 
  filter(plant_name_id %in% srli_raw$plant_name_id) 

srli_names_threat <- plants_full %>% 
  filter(plant_name_id %in% srli_threat$plant_name_id) 

prop_srli <- plantlist_dist_phylo_growth %>% 
  filter(plant_name_id %in% srli_names_raw$plant_name_id) %>% 
  group_by(LEVEL3_COD, LEVEL1_NAM) %>%  
  summarise(prop_org = n_distinct(plant_name_id) / nrow(srli_names_raw)) %>% 
  right_join(continent_names, by = c("LEVEL3_COD", "LEVEL1_NAM"))  %>% 
  mutate(prop_org = ifelse(is.na(prop_org), 0, prop_org))




# continental ------------------------------------------------------------------


species_samples_red_list <- readRDS("output/random_samples/redlist_10000_samples.rds")

redlist_null <- lapply(species_samples_red_list, subsampling.plants.null.cont, source = redlist_names_raw) %>% 
  rbindlist(idcol = "sample") %>% 
  mutate(dataset = "redlist")


species_samples_srli_list <- readRDS("output/random_samples/srli_10000_samples.rds")


srli_null <- lapply(species_samples_srli_list, subsampling.plants.null.cont, 
                    source = srli_names_raw) %>% 
  rbindlist(idcol = "sample") %>% 
  mutate(dataset = "srli")


# summary stats ----------------------------------------------------------------

redlist_null_stats <- redlist_null %>% 
  group_by(LEVEL1_NAM, dataset) %>% 
  summarise(   n = n(),
               mean = mean(prop_diff , na.rm = TRUE),
               sd = sd(prop_diff , na.rm = TRUE),
               se = sd / sqrt(n),
               ci_lower = mean - qt(0.975, df = n - 1) * se,
               ci_upper = mean + qt(0.975, df = n - 1) * se,
               q_lower = quantile(prop_diff, probs = 0.025), 
               q_upper = quantile(prop_diff, probs = 0.975)
  ) %>% 
  mutate(sig = !(q_lower <= 0 & q_upper >= 0))

srli_null_stats <- srli_null %>% 
  group_by(LEVEL1_NAM, dataset) %>% 
  summarise(n = n(),
            mean = mean(prop_diff , na.rm = TRUE),
            sd = sd(prop_diff , na.rm = TRUE),
            se = sd / sqrt(n),
            ci_lower = mean - qt(0.975, df = n - 1) * se,
            ci_upper = mean + qt(0.975, df = n - 1) * se,
            q_lower = quantile(prop_diff, probs = 0.025), 
            q_upper = quantile(prop_diff, probs = 0.975)
  ) %>% 
  mutate(sig = !(q_lower <= 0 & q_upper >= 0))


null_stats <- rbind(redlist_null_stats, srli_null_stats) 

null_stats_w <- null_stats %>% 
  dplyr::select(LEVEL1_NAM, mean,q_lower, q_upper, dataset) %>% 
  pivot_wider(
    names_from = dataset,
    values_from = c(mean,q_lower, q_upper)
  ) %>% 
  summarise(mean_diff =mean_srli- mean_redlist, 
            low_diff = q_lower_srli- q_lower_redlist,
            high_diff = q_upper_srli- q_upper_redlist,
            abs_diff = abs(mean_srli) - abs(mean_redlist)) %>% 
  mutate(sig = !(low_diff <= 0 & high_diff >= 0))


# calculate sample diff --------------------------------------------------------

# proportion of threatened species corrected to expected proportions under
# random sampling by diving proportion observed / proportion expected

nulls <- rbind(redlist_null, srli_null) %>% 
  mutate(prop_cor = n_prop / prop_org)

# sp_red <- dist_native_calc %>% 
#   filter(plant_name_id %in% redlist_names_raw$plant_name_id) %>% 
#   group_by(LEVEL1_NAM) %>%  
#   summarise(n_org = n_distinct(plant_name_id)) 
# 
# prop_red_threat <- dist_native_calc %>% 
#   filter(plant_name_id %in% redlist_names_threat$plant_name_id) %>% 
#   group_by(LEVEL1_NAM) %>%  
#   summarise(n_threat = n_distinct(plant_name_id), 
#             dataset = "redlist") %>% 
#   left_join(sp_red, by = "LEVEL1_NAM") %>% 
#   mutate(prop_threat =  n_threat /n_org, 
#          prop_org_tot = n_org / nrow(redlist_names_raw), 
#          prop_threat_tot = n_threat /nrow(redlist_names_raw))
# 
# sp_srli <- dist_native_calc %>% 
#   filter(plant_name_id %in% srli_names_raw$plant_name_id) %>% 
#   group_by(LEVEL1_NAM) %>%  
#   summarise(n_org = n_distinct(plant_name_id)) 
# 
# prop_threat_comb <- dist_native_calc %>% 
#   filter(plant_name_id %in% srli_names_threat$plant_name_id) %>% 
#   group_by(LEVEL1_NAM) %>%  
#   summarise(n_threat = n_distinct(plant_name_id), 
#             dataset = "srli") %>% 
#   left_join(sp_srli, by = "LEVEL1_NAM") %>% 
#   mutate(prop_threat =  n_threat / n_org, 
#          prop_org_tot = n_org / nrow(srli_names_raw), 
#          prop_threat_tot = n_threat /nrow(srli_names_raw) ) %>%  
#   rbind(prop_red_threat) 


plantlist_names <- plants_full %>%  
  dplyr::select(plant_name_id, taxon_rank, family, taxon_name, growth_form)

dist_native_calc_cont <- dist_native_calc %>% 
  group_by(LEVEL1_NAM, plant_name_id) %>% 
  summarise() %>% 
  left_join(plantlist_names, by = "plant_name_id") 


plantlist_dist_threat_red <- dist_native_calc_cont %>% 
  filter(plant_name_id %in% redlist_raw$plant_name_id) %>% 
  mutate(threat = case_when(
    plant_name_id %in% redlist_threat$plant_name_id     ~ "threatened",
    plant_name_id %in% redlist_dd$plant_name_id         ~ "data_deficient",
    plant_name_id %in% redlist_ext$plant_name_id         ~ "extinct",
   # plant_name_id %in% redlist_ntlc$plant_name_id         ~ "not_threatened",
    
    TRUE ~ "not_threatened"
  ))
table(plantlist_dist_threat_red$threat)
n_distinct(plantlist_dist_threat_red$plant_name_id)

plantlist_dist_threat_srli <- dist_native_calc_cont %>% 
  filter(plant_name_id %in% srli_raw$plant_name_id) %>% 
  left_join(plantlist_names, by = "plant_name_id") %>% 
  mutate(threat = case_when(
    plant_name_id %in% srli_threat$plant_name_id     ~ "threatened",
    plant_name_id %in% srli_dd$plant_name_id         ~ "data_deficient",
    plant_name_id %in% srli_ext$plant_name_id         ~ "extinct",
    
    TRUE ~ "not_threatened"
  ))
table(plantlist_dist_threat_srli$threat)
n_distinct(plantlist_dist_threat_srli$plant_name_id)

# calculations of threat counts 
sp_red <- dist_native_calc %>% 
  filter(plant_name_id %in% redlist_raw$plant_name_id) %>% 
  filter(!plant_name_id %in% redlist_ext$plant_name_id) %>% 
  group_by(LEVEL1_NAM) %>%  
  summarise(n_org = n_distinct(plant_name_id)) 

sp_red_dd <- dist_native_calc %>% 
  filter(plant_name_id %in% redlist_dd$plant_name_id) %>% 
  group_by(LEVEL1_NAM) %>%  
  summarise(n_dd = n_distinct(plant_name_id)) 


prop_red_threat <- dist_native_calc %>% 
  filter(plant_name_id %in% redlist_threat$plant_name_id) %>% 
  group_by(LEVEL1_NAM) %>%  
  summarise(n_threat = n_distinct(plant_name_id), 
            dataset = "redlist") %>% 
  left_join(sp_red, by = "LEVEL1_NAM") %>% 
  left_join(sp_red_dd, by = "LEVEL1_NAM") %>%
  replace(is.na(.), 0) %>% 
  mutate(
    prop_threat_lower =  n_threat / n_org, 
    prop_threat = n_threat / (n_org - n_dd), 
    prop_threat_upper =  (n_threat+n_dd) / n_org,
    prop_org_tot = n_org / nrow(srli_raw), 
    prop_threat_tot = n_threat /nrow(srli_raw)) 

sp_srli <- dist_native_calc %>% 
  filter(plant_name_id %in% srli_raw$plant_name_id) %>% 
  filter(!plant_name_id %in% srli_ext$plant_name_id) %>% 
  group_by(LEVEL1_NAM) %>%  
  summarise(n_org = n_distinct(plant_name_id)) 


sp_srli_dd <- dist_native_calc %>% 
  filter(plant_name_id %in% srli_dd$plant_name_id) %>% 
  group_by(LEVEL1_NAM) %>%  
  summarise(n_dd = n_distinct(plant_name_id)) 

prop_threat_comb <- dist_native_calc %>% 
  filter(plant_name_id %in% srli_threat$plant_name_id) %>% 
  group_by(LEVEL1_NAM) %>%  
  summarise(n_threat = n_distinct(plant_name_id), 
            dataset = "srli") %>% 
  left_join(sp_srli, by = "LEVEL1_NAM") %>% 
  left_join(sp_srli_dd, by = "LEVEL1_NAM") %>% 
  replace(is.na(.), 0) %>% 
  mutate(
    prop_threat_lower =  n_threat / n_org, 
    prop_threat = n_threat / (n_org - n_dd), 
    prop_threat_upper =  (n_threat+n_dd) / n_org,
    prop_org_tot = n_org / nrow(srli_raw), 
    prop_threat_tot = n_threat /nrow(srli_raw))  %>% 
  rbind(prop_red_threat)


# calculate expected threat ----------------------------------------------------
# correct threatened proportions
threat_null <- prop_threat_comb %>% 
  right_join(nulls, by = c("LEVEL1_NAM", "dataset")) %>% 
  mutate(prop_threat_cor = prop_threat * prop_cor)

threat_stats <- nulls %>%  
  left_join(prop_threat_comb %>%  
              dplyr::select(LEVEL1_NAM, prop_threat_est = prop_threat,
                            prop_threat_lower, prop_threat_upper, n_dd, dataset), 
            by = c("LEVEL1_NAM", "dataset")) %>% 
  mutate(prop_threat_upper_cor = prop_threat_upper * prop_cor, 
         prop_threat_est_cor = prop_threat_est * prop_cor, 
         prop_threat_lower_cor = prop_threat_lower * prop_cor) %>% 
  group_by(LEVEL1_NAM, dataset) %>% 
  summarise(
    mean = mean(prop_threat_est_cor), 
    
    q_lower = quantile(prop_threat_lower_cor, probs = 0.025), 
    q_upper = quantile(prop_threat_upper_cor, probs = 0.975), 
    dataset = unique(dataset))

threat_null_stats <- threat_null %>%
  group_by(LEVEL1_NAM, dataset) %>%
  summarise(n = n(),
            mean = mean(prop_threat_cor , na.rm = TRUE),
            sd = sd(prop_threat_cor , na.rm = TRUE),
            se = sd / sqrt(n),
            ci_lower = mean - qt(0.975, df = n - 1) * se,
            ci_upper = mean + qt(0.975, df = n - 1) * se,
            q_lower = quantile(prop_threat_cor, probs = 0.025),
            q_upper = quantile(prop_threat_cor, probs = 0.975)
  )




null_stats_w <- threat_stats %>% 
  dplyr::select(LEVEL1_NAM, mean,q_lower, q_upper, dataset) %>% 
  pivot_wider(
    names_from = dataset,
    values_from = c(mean,q_lower, q_upper)
  ) %>% 
  summarise(mean_diff =mean_redlist- mean_srli, 
            low_diff = q_lower_redlist- q_lower_srli,
            high_diff = q_upper_redlist- q_upper_srli,
            
            min_possible_gap = q_lower_redlist  - q_upper_srli, 
            max_possible_gap = q_upper_redlist - q_lower_srli,
            
            abs_diff = abs(mean_redlist) - abs(mean_srli)) %>% 
  mutate(sig = !(low_diff <= 0 & high_diff >= 0)) %>% 
  mutate(across(where(is.numeric), ~ round(.x * 100, 2))) 

# absolute difference for paper

red_all_threat <- plantlist_dist_threat_red %>% 
  group_by(plant_name_id, threat) %>% 
  summarise()

overall <- as.data.frame(table(red_all_threat$threat))
red_overall_tib <- tibble(
  n_org = sum(overall[overall$Var1 == "threatened",]$Freq + overall[overall$Var1 == "not_threatened",]$Freq + overall[overall$Var1 == "data_deficient",]$Freq), 
  n_threat = overall[overall$Var1 == "threatened",]$Freq, 
  n_dd = overall[overall$Var1 == "data_deficient",]$Freq,
  dataset="redlist"
  
)

prop_red_threat_all <- red_overall_tib %>% 
  mutate(
    prop_threat_lower =  n_threat / n_org, 
    prop_threat = n_threat / (n_org - n_dd), 
    prop_threat_upper =  (n_threat+n_dd) / n_org,
    prop_org_tot = n_org / nrow(srli_raw), 
    prop_threat_tot = n_threat /nrow(srli_raw),
    dataset="redlist"
  ) 


srli_all_threat <- plantlist_dist_threat_srli %>% 
  group_by(plant_name_id, threat) %>% 
  summarise()

overall <- as.data.frame(table(srli_all_threat$threat))
srli_overall_tib <- tibble(
  n_org = sum(overall[overall$Var1 == "threatened",]$Freq + overall[overall$Var1 == "not_threatened",]$Freq + overall[overall$Var1 == "data_deficient",]$Fre), 
  n_threat = overall[overall$Var1 == "threatened",]$Freq, 
  n_dd = overall[overall$Var1 == "data_deficient",]$Freq,
  dataset="srli"
  
)

prop_srli_threat_all <- srli_overall_tib %>% 
  mutate(
    prop_threat_lower =  n_threat / n_org, 
    prop_threat = n_threat / (n_org - n_dd), 
    prop_threat_upper =  (n_threat+n_dd) / n_org,
    prop_org_tot = n_org / nrow(srli_raw), 
    prop_threat_tot = n_threat /nrow(srli_raw), 
    dataset="srli") %>%  rbind(prop_red_threat_all)


stats_w <-  prop_srli_threat_all %>% 
  dplyr::select(prop_threat_lower, prop_threat,prop_threat_upper, dataset) %>% 
  pivot_wider(
    names_from = dataset,
    values_from = c(prop_threat_lower,prop_threat, prop_threat_upper)
  ) %>% 
  summarise(mean_diff = prop_threat_srli - prop_threat_redlist , 
            low_diff = prop_threat_lower_srli- prop_threat_lower_redlist,
            high_diff = prop_threat_upper_srli- prop_threat_upper_redlist,

            min_possible_gap = prop_threat_upper_srli  - prop_threat_lower_redlist, 
            max_possible_gap = prop_threat_lower_srli - prop_threat_upper_redlist ) %>% 
  mutate(sig = !(low_diff <= 0 & high_diff >= 0)) %>% 
  mutate(across(where(is.numeric), ~ round(.x * 100, 2)))



fwrite(nulls, "red_srli_null_prop_09_2026.txt")
fwrite(null_stats, "red_srli_null_prop_stats_09_2026.txt")
fwrite(threat_null, "red_srli_null_prop_threat_09_2026.txt")
fwrite(threat_null_stats, "red_srli_null_prop_stats_threat_09_2026.txt")
fwrite(threat_stats, "red_srli_threat_stats_cor_fact_09_2026.txt")


