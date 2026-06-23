library(tidyverse)
library(data.table)
library(sf)
library(parallel)
library(zoomerjoin)
library(rWCVP)
requireNamespace("rWCVPdata")

redlist_full <- read.csv("data/red/redlist_17_10_2025/taxonomy.csv") %>% 
  mutate(authority = str_replace(authority, "amp;", ""))


redlist_index <- read.csv("srli_redlist_match_04_26.csv") %>% 
  left_join(redlist_full, by = "scientificName")


redlist_analysis <- redlist_full %>% 
  dplyr::select(scientificName, authority) %>% 
  rbind(redlist_index %>%   
          dplyr::select(scientificName, authority)) %>% 
  mutate_all(na_if,"") %>% 
  filter(!duplicated(scientificName)) %>% 
  mutate(name_authority_org = ifelse(!is.na(authority), paste(scientificName, authority), scientificName))

length(unique(redlist_analysis$scientificName))


wcvp_raw <- fread("data/wcvp/wcvp_names_032023.csv", header = T)
dist_raw <- fread("data/wcvp/wcvp_distribution_032023.csv", header = T) 



plants_full_all <- fread("data/wcvp_accepted_merged.txt")
dist_native_all <- fread("data/dist_native.txt")


matched_names <- wcvp_match_names(names_df = redlist_analysis, wcvp_names =  wcvp_raw, 
                                  name_col = "scientificName",     
                                  author_col="authority",progress_bar = TRUE) 

names(matched_names) <- str_replace(names(matched_names), "plant_name_", "wcvp_")

xx <- matched_names %>% 
  group_by(scientificName) %>% 
  filter(n() > 1) %>% 
  slice_min(wcvp_author_edit_distance)


fuzzy_matching <- matched_names %>% 
  filter(!is.na(match_type)) %>% 
  filter(match_similarity > 0.8) %>% 
  group_by(scientificName) %>% 
  slice_max(match_similarity, with_ties =  T) %>% 
  slice_min(wcvp_author_edit_distance, with_ties = F)

xx <- fuzzy_matching %>% 
  filter(!match_similarity == 1)
table(xx$match_type)

length(unique(fuzzy_matching$wcvp_accepted_id))


wcvp_data_accepted <- fuzzy_matching %>% 
  filter(wcvp_rank == "Species") %>% 
  filter(wcvp_status == "Accepted") %>%  
  dplyr::select(plant_name_id = wcvp_id, taxon_name = wcvp_name, scientificName)

wcvp_data <- fuzzy_matching %>% 
  filter(!wcvp_status == "Accepted") %>% 
  filter(!scientificName %in% wcvp_data_accepted$scientificName)
  
names_no_acc_sp <- fuzzy_matching %>%  
    filter(!wcvp_status == "Illegitimate") %>% 
    filter(!scientificName %in% wcvp_data_accepted$scientificName) %>% 
    dplyr::select(plant_name_id_old = wcvp_id,taxon_name_old = wcvp_name, accepted_plant_name_id = wcvp_accepted_id, scientificName)
  
  accepted_ids_match <- wcvp_raw %>%
    #filter(taxon_status == "Accepted") %>% 
    filter(plant_name_id %in% names_no_acc_sp$accepted_plant_name_id) 
  
  accepted_sp <- accepted_ids_match %>%
    filter(taxon_status == "Accepted" & taxon_rank == "Species") %>% 
    filter(plant_name_id %in% names_no_acc_sp$accepted_plant_name_id) %>% 
    dplyr::select(taxon_name, plant_name_id) %>% 
    left_join(names_no_acc_sp, by = c("plant_name_id" = "accepted_plant_name_id"))
  

  not_sp <- accepted_ids_match %>% 
    filter(!taxon_rank == "Species") %>% 
    dplyr::select(plant_name_id_old_nosp = plant_name_id, parent_plant_name_id, 
                  taxon_name_old_nosp = taxon_name)
  
  accepted_sp_from_nosp <-  wcvp_raw %>%
    #filter(taxon_status == "Accepted") %>% 
    filter(plant_name_id %in% not_sp$parent_plant_name_id) %>% 
    dplyr::select(taxon_name, plant_name_id) %>% 
    left_join(not_sp, by = c("plant_name_id" = "parent_plant_name_id")) %>% 
    left_join(names_no_acc_sp, by = c("plant_name_id_old_nosp" = "accepted_plant_name_id"))
  
  
  ill_names <- fuzzy_matching %>%  
    filter(wcvp_status == "Illegitimate") %>% 
    ungroup() %>% 
    mutate(
      # taxon_name_old_nosp = case_when(!wcvp_rank == "Species" ~ wcvp_name,
      #                                TRUE ~ NA),
      taxon_name_old = case_when(wcvp_rank == "Species" ~ wcvp_name, 
                                 TRUE ~ NA), 
      # plant_name_id_old_nosp = case_when(!wcvp_rank == "Species" ~ wcvp_id,
      #                                TRUE ~ NA),
      plant_name_id_old = case_when(wcvp_rank == "Species" ~ wcvp_id, 
                                 TRUE ~ NA)
      
      ) %>% 
    dplyr::select(plant_name_id_old,taxon_name_old,
                  accepted_plant_name_id = wcvp_accepted_id, scientificName)
  
  accepted_legal_spid <- wcvp_raw %>%
    filter(plant_name_id %in% ill_names$accepted_plant_name_id) %>%  
    mutate(
      plant_name_id_new = case_when(taxon_rank == "Species" ~ accepted_plant_name_id, 
                                         taxon_rank != "Species" ~ parent_plant_name_id, 
                                    
           TRUE ~ NA), 
      taxon_name_old_nosp = case_when(!taxon_rank == "Species" ~ taxon_name,
                                      TRUE ~ NA),
      plant_name_id_old_nosp = case_when(!taxon_rank == "Species" ~ plant_name_id,
                                         TRUE ~ NA)) %>% 
    dplyr::select(taxon_name_old_nosp, plant_name_id_new,plant_name_id_old_nosp, accepted_plant_name_id) %>% 
    left_join(ill_names, by = "accepted_plant_name_id") %>% 
    mutate(plant_name_id = plant_name_id_new, 
           plant_name_id_old_nosp = case_when(plant_name_id_new != accepted_plant_name_id ~ accepted_plant_name_id, 
                                              TRUE ~ NA)) %>% 
    filter(!is.na(plant_name_id)) %>% 
    dplyr::select(-c(plant_name_id_new,accepted_plant_name_id))
    
    
    
  accepted_legal_sp <- wcvp_raw %>%
    filter(plant_name_id %in% accepted_legal_spid$plant_name_id)  %>% 
    dplyr::select(taxon_name, plant_name_id) %>% 
    right_join(accepted_legal_spid, by = "plant_name_id")
  
  
full_sp <- bind_rows(accepted_sp, accepted_sp_from_nosp, wcvp_data_accepted, accepted_legal_sp) 


length(unique(full_sp$taxon_name))



redlist_full_ass <- fread(("data/red/redlist_17_10_2025/assessments.csv")) %>% 
  mutate(year = as.numeric(substr(assessmentDate, 1, 4)))  %>% 
  mutate(redlistCategory = case_when(
    redlistCategory == "Lower Risk/conservation dependent" ~ "Least Concern", 
    redlistCategory ==  "Lower Risk/least concern" ~ "Least Concern", 
    redlistCategory ==  "Lower Risk/near threatened" ~ "Near Threatened", 
    T ~ redlistCategory))

redlist <- full_sp %>% 
  left_join(redlist_full_ass, by = "scientificName") %>% 
  filter(!is.na(redlistCategory))



levels_order <- c(
  "Least Concern",
  "Near Threatened",
  "Vulnerable",
  "Endangered",
  "Critically Endangered",
  "Extinct in the Wild",
  "Extinct", 
  "Data Deficient"
)



all_dup <- datawizard::data_duplicated(redlist, select = "taxon_name")


redlist_dups1 <- all_dup %>% 
  filter(taxon_name == scientificName) %>%
  ungroup() %>% 
  dplyr::select(-c(Row, count_na)) 

redlist_dups <- all_dup %>% 
  group_by(taxon_name) %>%
  filter(!any(taxon_name == scientificName, na.rm = TRUE)) %>%
  ungroup() %>%  
  mutate(status_ordered = factor(redlistCategory, levels = levels_order, ordered = TRUE))  %>% 
  group_by(taxon_name) %>%
  slice_min(status_ordered, with_ties = F) %>% 
  ungroup() %>%  
  dplyr::select(-c(Row, status_ordered, count_na)) %>% 
  rbind(redlist_dups1)

length(unique(redlist_dups$taxon_name)) - length(unique(all_dup$taxon_name))


redlist_full <- redlist %>% 
  filter(!taxon_name %in% redlist_dups$taxon_name) %>% 
  rbind(redlist_dups) %>% 
  filter(plant_name_id %in% plants_full_all$plant_name_id)

xx <- datawizard::data_duplicated(redlist_full, select = "plant_name_id")

redlist_full_adj <-  redlist_full 

table(redlist_full_adj$redlistCategory)


redlist %>% 
  filter(!taxon_name %in% redlist_full_adj$taxon_name) %>% 
  count()


write.csv(redlist_full_adj, "data/red/cleaned_10_2025/redlist_data_04_2026.csv")


# srli -------------------------------------------------------------------------

length(unique(redlist_index$scientificName))

srli_all <- redlist_index %>% 
  left_join(full_sp, by = "scientificName") %>% 
  filter(!is.na(plant_name_id))

srli_no_cat <- srli_all %>% 
  filter(is.na(redlistCategory)) %>% 
  left_join(redlist_full_ass %>%  
              dplyr::select(red_cat = redlistCategory, 
                                                red_year = year, 
                                                red_year_pub = yearPublished, 
                                                scientificName), by = "scientificName") %>%  
  filter(red_cat %in% unique(srli_all$redlistCategory)) %>% 
  mutate(redlistCategory = red_cat, 
         yearPublished = red_year_pub, 
         year = red_year
  )  %>% 
  dplyr::select(any_of(names(srli_all)))

sum(table(srli_no_cat$redlistCategory))

srli <- srli_all %>% 
  filter(!scientificName %in% srli_no_cat$scientificName) %>% 
  rbind(srli_no_cat)


all_dup_srli <- datawizard::data_duplicated(srli, select = "taxon_name")


srli_dups1 <- all_dup_srli %>% 
  filter(taxon_name == scientificName) %>%
  ungroup() %>% 
  dplyr::select(-c(Row, count_na)) 

srli_dups <- all_dup_srli %>% 
  group_by(taxon_name) %>%
  filter(!any(taxon_name == scientificName, na.rm = TRUE)) %>%
  ungroup() %>%  
  mutate(status_ordered = factor(redlistCategory, levels = levels_order, ordered = TRUE))  %>% 
  group_by(taxon_name) %>%
  slice_min(status_ordered, with_ties = F) %>% 
  ungroup() %>%  
  dplyr::select(-c(Row, status_ordered, count_na)) %>% 
  rbind(srli_dups1)


srli_full <- srli %>% 
  filter(!taxon_name %in% srli_dups$taxon_name) %>% 
  #rbind(srli_dups) %>% 
  filter(plant_name_id %in% plants_full_all$plant_name_id)  %>% 
  filter(!is.na(redlistCategory))

n_distinct(srli_full$plant_name_id) - nrow(srli_full)
table(srli_full$redlistCategory)

write.csv(srli_full, "data/red/cleaned_10_2025/srli_data_04_2026.csv")

# 
# redlist_in_srli <- redlist_full_adj  %>% 
#   filter(plant_name_id %in% srli_full$plant_name_id) %>% 
#   filter(year > 2015)
# 
# table(redlist_in_srli$redlistCategory) / sum(table(redlist_in_srli$redlistCategory))
# 
# 
# srli_in_red <- srli_full  %>% 
#   filter(plant_name_id %in% redlist_in_srli$plant_name_id) %>% 
#   filter(year < 2015)
# 
# table(srli_in_red$redlistCategory) / sum(table(srli_in_red$redlistCategory))
# 
# redlist_match_srli <- redlist_full_adj  %>% 
#   filter(plant_name_id %in% srli_in_red$plant_name_id)
# 
# table(redlist_match_srli$redlistCategory) / sum(table(redlist_match_srli$redlistCategory))
