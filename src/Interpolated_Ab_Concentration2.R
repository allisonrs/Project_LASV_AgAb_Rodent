##################### ELISA IgG Cutoffs ##################################
### 4PL curve of mean ODs from Prism --> starting data is interpolated files ("data/ELISA_interpolation_csv/"). 
###  The standard curve (3x dilution) per plate is interpolated from 1000U/mL (so staring conc is log(1000)=3) 

library(tidyverse)

# Compile all Interpolated IgG ELISA data
Interpolated <- list.files(path = "data/ELISA_interpolation_csv/", pattern = "Interpolation.*\\.csv", full.names = TRUE) %>%
  map_df(~read_csv(., col_types = "c", id = "File")) %>%
  mutate(Species_ab = if_else(str_detect(File, "Mouse"), "Mouse", "Rat")) %>%
  mutate(Protein = if_else(str_detect(File, "GP"), "GP", "NP")) %>% 
  select(Animal_ID, Species_ab, Protein, U_mL) 


# IgG of animals only -- removing blanks and negative controls from each plate
Interpolated_Animals <- Interpolated %>% 
  filter(!str_detect(Animal_ID, "[:upper:]")) %>%
  mutate(Animal_ID = as.numeric(unlist(.$Animal_ID))) #%>%
  #write_csv(file = "data/LASV_IgG_Animals_Interpolated.csv") 


# Negative Cutoffs as determined by known negative serum controls (ultimately not used for cutoff)
Neg <- list.files(path = "data/ELISA_interpolation_csv/Negatives/", pattern = "Interpolation.*\\.csv", full.names = TRUE) %>%
  map_df(~read_csv(., id = "File")) %>% 
  mutate(Species_ab = if_else(str_detect(File, "Rat"), "Rat", "Mouse")) %>%
  mutate(Protein = if_else(str_detect(File, "GP"), "GP", "NP")) %>%
  select(Species_ab, Protein, U_mL) %>% 
  group_by(Species_ab, Protein) %>% 
  summarize(ave_neg_U_mL = mean(U_mL, na.rm = TRUE), "2x_ave_neg_U_mL" = 2*(mean(U_mL, na.rm = TRUE)), n=sum(!is.na(U_mL))) 

# Negative cutoffs (2x OD of Negative Control or Blank of experimental plate, whichever is higher) across all experimental plates 
Negatives <- Interpolated %>%
  group_by(Species_ab, Protein) %>% 
  filter(str_detect(Animal_ID, "2x.*")) %>%
  summarize(ave_neg_U_mL = mean(U_mL, na.rm = TRUE), n=sum(!is.na(U_mL))) 


# Determine which samples are positive based on Negative cutoffs calculated directly above
Interpolated_Animals  %>%
  left_join(Negatives, by = c("Species_ab", "Protein"))%>% 
  select(1:5) %>%
  mutate(Positive = (U_mL >= ave_neg_U_mL)) %>%
  replace_na(list(Positive = FALSE)) #%>% #some of the tested animals had such low ODs that it could not be interpolated
  #write.csv(file = "data/LASV_IgG_Animals_PosNeg.csv")


# Append ELISA results to rest of data to create main analysis 
Animal_data <- read_csv(file = "data/LASV_Rodent_data.csv")
Interpolated_Animals_Results <- read.csv(file = "data/LASV_IgG_Animals_PosNeg.csv")
  
Interpolated_Animals_Results %>%
  select(Animal_ID, Protein, U_mL, Positive) %>%
  pivot_wider(id_cols = "Animal_ID", names_from = Protein, values_from = c(U_mL, Positive)) %>%
  unnest(cols = c(U_mL_GP, U_mL_NP, Positive_GP, Positive_NP)) %>%
  full_join(Animal_data, by = "Animal_ID") %>% 
  arrange(Animal_ID) %>% 
  mutate(IgG_Positive = case_when(Positive_GP == TRUE | Positive_NP == TRUE ~ TRUE,
                                  Positive_GP == FALSE & Positive_NP == FALSE ~ FALSE,
                                  TRUE ~ NA)) %>% 
  relocate(2:5, .after = RDT_Score) %>% 
  relocate(Positive_GP, .after = U_mL_GP) #%>% 
  #write_csv(file = "data/LASV_Analysis.csv")
