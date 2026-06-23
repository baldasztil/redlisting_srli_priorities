library(data.table)
library(tidyverse)
library(performance)
library(sf)


# threat correction ------------------------------------------------------------

nulls <- fread("output/redlist_srli/plotting_data/red_srli_null_prop_threat_04_2026.txt")

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


fwrite(threat_stats, "red_srli_threat_stats_cor_fact_04_2026.txt")


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

