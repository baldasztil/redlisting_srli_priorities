library(tidyverse)
library(data.table)
library(sf)
library(parallel)
library(zoomerjoin)

redlist_analysis <- read.csv("data/red/redlist_17_10_2025/assessments.csv") %>% 
  filter(!redlistCategory %in% c(
    "Extinct in the Wild",
    "Extinct",
    "Not Evaluated")) %>% 
  mutate(year =  as.numeric(substr(assessmentDate, 1, 4)))

redlist_syns <- read.csv("data/red/redlist_17_10_2025/synonyms.csv") %>% 
  mutate(sp1 = paste(genusName, speciesName))  %>% 
  filter(!sp1 %in% redlist_analysis$scientificName)

pterido <- fread("PteridophyteRatingsForLudwig.csv") %>% 
  mutate(redlistCategory = case_when(rl_category == "LC" ~ "Least Concern", 
                                     rl_category == "EN" ~ "Endangered", 
                                     rl_category == "VU" ~ "Vulnerable", 
                                     rl_category == "CR" ~ "Critically Endangered", 
                                     rl_category == "DD" ~ "Data Deficient", 
                                     rl_category == "NT" ~ "Near Threatened")) %>% 
  mutate(id = "Pteridophytes", 
    year = NA, 
    yearPublished = NA) %>% 
  dplyr::select(sp1 = Rank, scientificName = Rank , redlistCategory, year, yearPublished, id) 

levels_order <- c(
  "Least Concern",
  "Near Threatened",
  "Vulnerable",
  "Endangered",
  "Critically Endangered",
  "Data Deficient"
)

redlist_historical <- read.csv("historical_assessments.csv") %>% 
  filter(!redlistCategory %in% c(
    "Extinct in the Wild",
    "Extinct",
    "Not Evaluated")) %>% 
  filter(between(year, 2001, 2015)) %>% 
  group_by(scientificName, redlistCategory, year, yearPublished) %>% 
  summarise() %>% 
  group_by(scientificName) %>% 
  slice_max(year, with_ties = T) %>% 
  mutate(status_ordered = factor(redlistCategory, levels = levels_order, ordered = TRUE))  %>% 
  group_by(scientificName) %>%
  slice_min(status_ordered, with_ties = F) 


cor_tax_srli <- fread("output/redlist_srli/srli/srli_noinfo_taxonomy_corr.csv",  header = T) %>% 
  dplyr::select(family,  genus,   species, sp1 = taxon_name,  id)  
#Get the path of filenames

list_data <- list.files("C:/Users/s2117440/OneDrive - Royal Botanic Garden Edinburgh/A_PhD/Chapters/Chapter1_UNEP_KEW/Projects/SubsamplingPlantBiodiversity/data/red/srli/csv/", full.names = TRUE) 
  

#Read them in a list
redlist_index <- lapply(list_data, fread) %>% 
  rbindlist()  %>%
  mutate(species = ifelse(is.na(species), "missing", species)) %>% 
  filter(!paste(genus, species)  %in%  
           paste(cor_tax_srli$genus, cor_tax_srli$species)) %>% 
  rbind(cor_tax_srli) %>% 
  filter(!duplicated(sp1)) %>% 
  filter(!sp1 %in% pterido$sp1)



wcvp_raw <- fread("data/wcvp/wcvp_names_032023.csv", header = T)
dist_raw <- fread("data/wcvp/wcvp_distribution_032023.csv", header = T) 

plants_full_all <- fread("data/wcvp_accepted_merged.txt")
dist_native_all <- fread("data/dist_native.txt")


# match syns 
matched_syns <- redlist_index %>% 
  filter(sp1 %in% redlist_syns$sp1) %>%  
  left_join(redlist_syns,  by = "sp1") %>% 
  dplyr::select(sp1,scientificName) %>% 
  left_join(redlist_analysis,  by = c("scientificName"))  %>% 
  dplyr::select(sp1,scientificName, redlistCategory, year, yearPublished)  
  

# filter all matches  
matched_acc_sp <- redlist_index %>% 
  filter(sp1 %in% redlist_analysis$scientificName) %>%  
  left_join(redlist_analysis,  by = c("sp1" = "scientificName"), keep = T) %>% 
  dplyr::select(sp1,scientificName, redlistCategory, year, yearPublished) 



matched_names_acc_syn <- rbind(matched_syns, matched_acc_sp) #%>% 
 # mutate(scientificName = ifelse(str_count(scientificName, " ") >= 2, str_extract(scientificName, "^\\S+\\s+\\S+"), scientificName)) 

# matching historical
matched_historical <- redlist_historical %>% 
  filter(scientificName %in% matched_names_acc_syn$scientificName) %>%  
  dplyr::select(scientificName, redlistCategory, year, yearPublished) %>% 
  left_join(matched_names_acc_syn %>% dplyr::select(scientificName, sp1), by = "scientificName") %>% 
  rbind(matched_names_acc_syn)

# matching historical missing
hist_missing_added <- redlist_historical %>% 
  filter(scientificName %in% redlist_index$sp1) %>% 
  filter(!scientificName %in% matched_historical$sp1) %>%  
  dplyr::select(scientificName, redlistCategory, year, yearPublished) %>% 
  left_join(matched_names_acc_syn %>% dplyr::select(scientificName, sp1), 
            by = "scientificName") %>% 
  rbind(matched_historical)  %>%  
  mutate(sp1 = ifelse(is.na(sp1), scientificName, sp1))  %>% 
  filter(between(year, 2001, 2015)) 


xx <- hist_missing_added %>% 
  group_by(scientificName) %>% 
  filter(n() > 1)
length(unique(hist_missing_added$scientificName))

matched_sp <- hist_missing_added$sp1



# filter time 
srli_assessments_red <- hist_missing_added %>% 
  #filter(between(year, 2001, 2015)) %>% 
  group_by(scientificName) %>% 
  slice_max(year, n =1, with_ties = F)

table(srli_assessments_red$redlistCategory)

# insert Pteridophytes ---------------------------- 

srli_assessments <- srli_assessments_red %>% 
  left_join(redlist_index %>% dplyr::select(sp1, id), by= "sp1") %>% 
  # filter(!redlistCategory == "Not Evaluated") %>% 
  # filter(!scientificName %in% pterido$scientificName) %>% 
  rbind(pterido) %>% 
  as.data.frame() 

missing <- redlist_index %>%   
  filter(!sp1 %in% srli_assessments$sp1) 

srli_full <- srli_assessments %>% 
  rbind(missing %>%
          mutate(
            scientificName = sp1, 
            year = NA, 
            yearPublished = NA, 
            redlistCategory = NA, 
            id = id) %>% 
          dplyr::select(scientificName, year, yearPublished, redlistCategory, sp1, id))


write.csv(srli_full, "srli_redlist_match_04_26.csv")

