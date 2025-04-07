library(tidyverse)
library(lubridate)
library(sf)
library(janitor)
library(gtsummary)
library(gt)
library(rstatix)
library(ggpubr)
library(ggprism)
library(epiR)
library(statpsych)
library(exact2x2)
library(boot)

### Prepare Data ##############################################################################

original_data <- read_csv("data/LASV_Analysis.csv") %>%
                  mutate(Date=dmy(.$Date))

analysis <- original_data %>%
  filter(!is.na(Ag_RDT), !is.na(IgG_Positive)) ## excludes specimens not tested

trap_nights <- read_csv("data/trap_nights_total_clean.csv") %>%
                mutate(date_set=mdy(.$date_set), 
                       date_recovered1=mdy(.$date_recovered1), 
                       date_recovered2=mdy(.$date_recovered2))

### Shape Files ##########

## administrative boundaries: https://data.humdata.org/dataset/cod-ab-sle

# Shape File of SL Provinces
sle_adm1 <- read_sf("data/SL_Map_Data/sle_admbnda_adm1_1m_gov_ocha_20161017.shp") %>%   
  clean_names()

# Shape File of SL Districts
sle_adm2 <- read_sf("data/SL_Map_Data/sle_admbnda_adm2_1m_gov_ocha_20161017.shp") %>% 
  clean_names() %>%
  filter(admin2name %in% "Kenema") ## filter to include Kenema district only

sle_adm3 <- read_sf("data/SL_Map_Data/sle_admbnda_adm3_1m_gov_ocha_20161017.shp") %>% 
  clean_names() %>%
  filter(admin2name %in% "Kenema") ## filter to include Kenema district only

## cities: 

# Shape file of SL cities (https://data.humdata.org/dataset/sierra-leone-settlements)

sle_pplp_nga <- read_sf("data/SL_Map_Data/sle_pplp_nga_20170523.shp") %>% 
  clean_names() %>%
  filter(feature_nam %in% "Freetown" | (feature_nam %in% "Kenema" & admin3name %in% "Kenema Town"))

kenema_town <- read_sf("data/SL_Map_Data/sle_pplp_nga_20170523.shp") %>% 
  clean_names() %>%
  filter(feature_nam %in% "Kenema" & admin3name %in% "Kenema Town")

### Tables and Figures in code (in order of appearance) ###############################

# Table A1: Trapping summary
# Table A2: Specimen counts by Cytochrome b PCR
# Table A3: Contingency Table of Ag/qRT-PCR status of serum of Ag+/- specimens
# Table A4: Sensitivity and specificity of serum Ag to qRT-PCR 
# Table A5: qRT-PCR status of select organs from select Ag+ specimens
# Figure 1: Map of trapping locations
# Table 1: Counts and row percentage of serology tests (Ag+, IgG+, IgG+ breakdown)
# Table 2: Contingency table of Ag/IgG +/- for Mastomys
# Table A6: Contingency table of Ag/IgG +/- for Rattus
# Figure 2: GP and NP IgG relative values by protein and genus
# Figure 3: Mean IgG relative value vs RDT Score for Mastomys, by protein
# Figure A1: Correlation between GP IgG vs NP IgG relative value in Mastomys



## Table A1: Trapping summary ##########################################################
# GPS_Coordinates of communities in study
gps_coordinates <- trap_nights %>% 
  group_by(community) %>%
  summarize(first(Lat_Convert), first(Long_Use)) %>%
  rename(lat = `first(Lat_Convert)`, lon = `first(Long_Use)`) 


# Trapping summary
trap_summary <-trap_nights %>%
  mutate(sum_trap_nights = case_when(is.na(date_recovered2) & is.na(date_recovered1) ~ 0,
                                     is.na(date_recovered2) ~ as.numeric(date_recovered1-date_set), 
                                     TRUE ~ as.numeric(date_recovered2-date_set)),
         sum_specimens = if_else(is.na(animal_id), 0, 1)) %>%
  group_by(community,section) %>%
  summarize(Duration = max(sum_trap_nights),
            N_homes=max(house_number),
            Total_trap_nights = sum(sum_trap_nights),
            Total_specimens = sum(sum_specimens)) %>%
  group_by(community) %>%
  summarize(duration = sum(Duration), 
            n_homes = sum(N_homes),
            total_trap_nights = sum(Total_trap_nights), 
            total_specimens = sum(Total_specimens))

# Counts of crocidura in community to exclude because off target species
crocidura <- original_data %>%
  mutate(Crocidura = if_else(Genus_og == 'Crocidura', 1, 0)) %>%
  group_by(Community) %>%
  count(Crocidura) %>% 
  filter(Crocidura ==1) %>%
  select(c(Community, n)) %>%
  rename(n_crocidura = 'n')


# Counts of non-crocidura specimens without ag or ab testing filtered from analysis
testing <- original_data %>%
  filter(Genus_og != "Crocidura") %>%
  mutate(Missing_tests = if_else(is.na(Ag_RDT) | is.na(IgG_Positive), 1, 0)) %>% 
  group_by(Community) %>%
  count(Missing_tests) %>%
  filter(Missing_tests == 1) %>%
  select(c(Community,n)) %>%
  rename(missing_tests = "n") 

## Putting it all together to create Table A1
Table_A1 <- left_join(gps_coordinates, trap_summary, by = "community") %>% 
  left_join(crocidura, by = c("community" = "Community")) %>% 
  mutate_if(is.numeric,coalesce,0) %>% 
  left_join(testing, by = c("community" = "Community")) %>% 
  mutate_if(is.numeric,coalesce,0) %>% 
  mutate(Total_included = total_specimens - (n_crocidura + missing_tests))

#write_csv(Table_A1, "Table_A1.csv")


### Table A2: Specimen counts by Cytochrome b PCR ####################################################

# Field-identified genus
collected <- analysis %>%
  group_by(Genus_og) %>%
  count() %>%
  arrange(desc(n)) %>%
  rename(Field_ID_Genus = Genus_og, Collected = n) 

# Field-identified genus ID of specimens that were tested by PCR
species_tested <- analysis %>%
  filter(!is.na(Cytb_species)) %>%
  group_by(Genus_og) %>%
  count() %>%
  arrange(desc(n)) %>%
  rename(Field_ID_Genus = Genus_og, Tested_by_PCR = n) 

# Species cytb test results
cytb <- analysis %>%
  group_by(Genus_og, Cytb_species) %>%
  count() %>%
  arrange(desc(n)) %>%
  rename(Field_ID_Genus = Genus_og) %>%
  filter(!is.na(Cytb_species))
  

## Putting it all together to create Table A2
Table_A2 <- left_join(collected, species_tested, by = "Field_ID_Genus") %>% 
  left_join(cytb, by = "Field_ID_Genus") %>%
  group_by(Field_ID_Genus) %>% 
  mutate(Collected = ifelse(row_number() == 1, Collected, NA),
         Tested_by_PCR = ifelse(row_number() == 1, Tested_by_PCR, NA)) %>% 
  ungroup()

write_csv(Table_S2, "Table_A2.csv")

### Table A3: Contingency Table of Ag/qRT-PCR status of serum of Ag+/- specimens ############################################################
# Serum qRT-PCR testing results (+/-)
analysis %>%
  filter(!is.na(qPCR_Positive)) %>%
  drop_na(Ag_RDT, Positive_AVL) %>%
  count(Positive_AVL)

# Frequency table of PCR/Ag results
freqtable_RDTPCR_serum <- analysis %>%
  filter(!is.na(qPCR_Positive)) %>%
  drop_na(Ag_RDT, Positive_AVL) %>% 
  mutate(Ag_RDT = if_else(Ag_RDT == TRUE, "Ag+", "Ag-"),
         qRTPCR = if_else(Positive_AVL == TRUE, "qRT-PCR+", "qRT-PCR-")) %>%
  mutate(Ag_RDT = factor(Ag_RDT, levels = c("Ag+", "Ag-")),
         qRTPCR = factor(qRTPCR, levels = c("qRT-PCR+", "qRT-PCR-"))) %>%
  xtabs(~Ag_RDT+qRTPCR, data = .) 

# Mcnemar_test
freqtable_RDTPCR_serum %>% mcnemar.exact()

## Output of Serum specificity/sensitivity and Menemar_test used to create Table A3.
## Code below contains Table A3 information w/o Odds ratio + 95% CI

freqtable_serum_RDTPCR_gt <- freqtable_RDTPCR_serum %>%
  as.data.frame.matrix() %>%
  mutate("Row Total" = rowSums(.)) %>%
  bind_rows(summarise(.,across(where(is.numeric), sum))) %>%
  rownames_to_column(var = " ") %>%
  assign_in(list(1, 3), "Column Total") %>%
  gt() %>%
  cols_label("qRT-PCR+" = md("**qRT-PCR +**"), "qRT-PCR-" = md("**qRT-PCR -**")) %>%
  tab_style(cell_text(weight = "bold"), location = cells_body(columns = " ", rows = 1:2)) %>%
  tab_source_note(source_note = "p = < 0.0001 by McNemar's Test") #%>% from above
#gtsave("freqtable_serum_RDTPCR_gt.png")

### Table A4: Sensitivity and specificity of serum Ag to qRT-PCR  #########################################################

# Serum specificity/sensitivity compared to qRT-PCR
serum_SS_RDTPCR_gt <- epi.tests(freqtable_RDTPCR_serum) %>% #From Table A3 above
  summary() %>%
  round(digits = 3) %>%
  rownames_to_column(var = "Characteristic") %>%
  as.data.frame.matrix() %>%
  .[3:4,] %>%
  mutate(Characteristic = if_else(Characteristic == "se", "Sensitivity", "Specificity")) %>%
  rename(" " = Characteristic) %>%
  gt() %>%
  cols_merge_range(col_begin = lower, col_end = upper, sep = "-") %>%
  fmt_percent(columns = c("est","lower", "upper"), decimals = 2) %>%
  cols_label(est = md("**est**"), lower = md("**95% CI**")) %>%
  tab_style(cell_text(weight = "bold"), location = cells_body(columns = " ", rows = 1:2)) #%>%
#gtsave("serum_SS_RDTPCR_gt.png")

### Table A5: qRT-PCR status of select organs from select Ag+ specimens ############################################################
# Filter for specimens with organs tested
PCR_tissues <- analysis %>%
  filter(!is.na(qPCR_Positive), Ag_RDT == T) %>% 
  mutate(AVL_only = if_else(!is.na(Positive_AVL) 
                            & is.na(Positive_Liver)
                            & is.na(Positive_Lung) 
                            & is.na(Positive_Kidney) 
                            &is.na(Positive_Spleen), 
                            "TRUE", "FALSE")) %>% 
  filter(AVL_only == F)

# At least one organ positive
PCR_tissues %>%
  group_by(qPCR_Positive) %>%
  count()

# Organ counts
PCR_tissues_counts <- PCR_tissues %>%
  rename_with(~gsub("Positive_", "", .x)) %>%
  pivot_longer(cols = c(AVL, Liver, Lung, Kidney, Spleen), names_to ="organ", values_to = "PCR_Positive_organ", values_drop_na = TRUE) %>% 
  filter(organ != "AVL") 

# Organ positivity
PCR_tissues_counts%>%
  group_by(organ) %>%
  count(PCR_Positive_organ) %>%
  pivot_wider(names_from = PCR_Positive_organ, values_from = n)

## Counts from above used to create Table A5

### Figure 1: Map of trapping locations  ###############

# Plot Prep
gps_sf<- st_as_sf(gps_coordinates, coords = c("lon", "lat"), crs = 4326)
st_crs(sle_adm3) <- 4326
st_crs(sle_adm2) <- 4326
intersected_shapes <- st_intersection(sle_adm3, sf_gps)

sle_adm_filtered <- sle_adm3 %>%
  filter(sle_adm3$admin3name %in% intersected_shapes$admin3name)

bbox <- st_bbox(sle_adm_filtered)


# Plot map of SL
SL <- ggplot() + 
  geom_sf(data = sle_adm1, fill = NA) + 
  geom_sf(data = sle_adm2, mapping = aes(fill = admin2name), fill = "grey") +
  geom_sf(data= sle_adm_filtered, mapping = aes(color = admin3name), color = "darkgrey") +
  geom_point(data = sle_pplp_nga, mapping = aes(geometry = geometry), 
             stat = "sf_coordinates", color = "black", size = 1, shape = if_else(sle_pplp_nga$feature_nam == "Freetown", 8, 18)) + #from coordinates in Table A1
  geom_sf_label(data = sle_pplp_nga,mapping = aes(label = feature_nam), size = 1.7, 
               nudge_x = if_else(sle_pplp_nga$feature_nam == "Freetown", 0.05, -0.25),
               nudge_y = if_else(sle_pplp_nga$feature_nam == "Freetown", -0.15, -0.05)) +
  geom_point(data = gps_coordinates, mapping = aes(x = lon, y = lat), shape = 20, size = 0.5) +
  geom_rect(aes(xmin = bbox["xmin"], xmax = bbox["xmax"], ymin = bbox["ymin"], ymax = bbox["ymax"]),
            color = "black", fill = NA, size = 0.5) +
  coord_sf() +
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'), 
        legend.key = element_rect(fill = 'white', color = 'white'), 
        legend.background = element_rect(fill = "white", color = "white"),
        axis.text.x = element_text(angle = 30))

# Plot the individual communities
zoomed <- ggplot() + 
  geom_sf(data = sle_adm_filtered, mapping = aes(fill = admin3name), fill = "lightgrey") +
  geom_sf(data = kenema_town, mapping = aes(color = feature_nam), color = "black", size = 2, shape = 18) +
  geom_sf_text(data = kenema_town,mapping = aes(label = feature_nam), size = 3,
               nudge_y = -0.04) +
  geom_point(data = gps_coordinates, mapping = aes(x = lon, y = lat), size = 1) +
  geom_text_repel(data = gps_coordinates, mapping = aes(x = lon, y = lat, label = community), 
                  max.overlaps = 20, min.segment.length = 0.4, size = 2, seed = 20) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'), 
        legend.key = element_rect(fill = 'white', color = 'white'), 
        legend.background = element_rect(fill = "white", color = "white"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

# Put the plots together
merge <- plot_grid(SL, zoomed, nrow = 1, rel_widths = c(2,1.5), rel_heights = c(2,0.5)) 
ggdraw(merge) +
  draw_line(x = c(0.46, 0.65), y = c(0.49, 0.55), arrow = arrow(length = unit(0.1, "inches")), color = "black") #+
  #ggsave("Figure1.eps", width = 2600, height = 1950, units = "px", dpi = 320) 
    
  
### Correct data frame for the one Malacomys later determined to be a Praomys ##########
  analysis[6,10] <- "Praomys" 


### Table 1: Counts and row percentage of serology tests (Ag+, IgG+, IgG+ breakdown) ########################################################################

# Function to calculate % positive for GP only or NP only
calc_pos_percent <- function(data, pos_col, neg_col) {
  pos_count <- sum(data[[pos_col]] == TRUE & data[[neg_col]] == FALSE)
  pos_percent <- round((pos_count / nrow(data)) * 100, 0)
  paste0(pos_count, " (", pos_percent, ")")
}

# All serology test results (n, %)
Table_1 <- analysis %>%
  rename(Genus = Genus_og) %>%
  group_by(Genus) %>%
  summarise(N = n(),
            Ag_n_percent = paste0(sum(Ag_RDT), " (", round((sum(Ag_RDT) / n()) * 100, 0), ")"),
            IgG_n_percent = paste0(sum(IgG_Positive), " (", round((sum(IgG_Positive) / n()) * 100, 0), ")"),
            GP_n_percent = calc_pos_percent(cur_data(), "Positive_GP", "Positive_NP"),
            NP_n_percent = calc_pos_percent(cur_data(), "Positive_NP", "Positive_GP"),
            GP_NP_n_percent = paste0(sum(Positive_GP == TRUE & Positive_NP == TRUE), " (", round((sum(Positive_GP == TRUE & Positive_NP == TRUE) / n()) * 100, 0), ")")) %>%
  arrange(desc(N))

#write.csv(Table_1, "Table_1.csv")


### Table 2 (Mastomys Contingency Table) and Table A6 (Rattus Contingency Table) ###############

# Mastomys contingency table
analysis %>%
  filter(Genus_og == "Mastomys") %>%
  drop_na(Ag_RDT, IgG_Positive) %>%
  mutate(Ag_RDT = if_else(Ag_RDT == TRUE, "Ag +", "Ag -"),
         IgG_Positive = if_else(IgG_Positive == TRUE, "IgG +", "IgG -")) %>%
  xtabs(~Ag_RDT+IgG_Positive, data = .) %>%
  addmargins()

# Rattus contingency table
analysis %>%
  filter(Genus_og == "Rattus") %>%
  drop_na(Ag_RDT, IgG_Positive) %>%
  mutate(Ag_RDT = if_else(Ag_RDT == TRUE, "Ag +", "Ag -"),
         IgG_Positive = if_else(IgG_Positive == TRUE, "IgG +", "IgG -")) %>%
  xtabs(~Ag_RDT+IgG_Positive, data = .) %>% #mcnemar.exact()   stop here for calculation
  #addmargins() 

# Contingency tables and p-values McNemar's 
analysis %>%
  filter(Genus_og == "Mastomys" | Genus_og == "Rattus") %>%
  mutate(Ag_RDT = if_else(Ag_RDT == TRUE, "Ag +", "Ag -"),
         IgG = if_else(IgG_Positive == TRUE, "IgG +", "IgG -")) %>%
  group_by(Genus_og) %>%
  select(Genus_og, Ag_RDT, IgG) %>%
  nest() %>%
  mutate(contingency = map(data, ~xtabs(~Ag_RDT+IgG, data = .)),
         x2 = map(contingency, mcnemar.exact), 
         glanced = map(x2, ~tibble::tibble(
           statistic = .$statistic,
           parameter = .$parameter,
           p.value = .$p.value,
           conf.int = .$conf.int,
           estimate = .$estimate,
           null.value = .$null.value,
           alternative = .$alternative,
           method = .$method,
           data.name = .$data.name
         ))) %>%
  unnest(c(glanced))

### Figure 2: GP and NP IgG relative values by protein and genus ########################################################

# Dataframes for analysis
analysis2 <- analysis %>% ## only IgG positive for at least one protein
  filter(IgG_Positive == "TRUE") %>%
  pivot_longer(cols = c("U_mL_GP", "U_mL_NP"), names_to = "Protein", values_to = "U_mL") %>% 
  mutate(Protein = if_else(Protein == "U_mL_GP", "GP", "NP")) %>%
  filter(Positive_GP == "TRUE" & Protein == "GP" | Positive_NP == "TRUE" & Protein == "NP") %>%
  mutate(Protein = fct_relevel(Protein, c("NP", "GP"))) %>%
  add_column(Outlier = FALSE) %>%
  mutate(Positive_GP = as.character(Positive_GP),
         Positive_NP = as.character(Positive_NP))

analysisGP <- analysis %>%
  pivot_longer(cols = c("U_mL_GP", "U_mL_NP"), names_to = "Protein", values_to = "U_mL") %>% 
  mutate(Protein = if_else(Protein == "U_mL_GP", "GP", "NP")) %>% 
  replace_na(list(U_mL = 0)) %>%
  filter(Genus_og == "Mastomys" | Genus_og == "Rattus", Protein =="GP") %>%
  add_column("Outlier" = FALSE) %>%
  mutate(Positive_GP = as.character(Positive_GP))
  

analysisNP <- analysis %>%
  pivot_longer(cols = c("U_mL_GP", "U_mL_NP"), names_to = "Protein", values_to = "U_mL") %>% 
  mutate(Protein = if_else(Protein == "U_mL_GP", "GP", "NP")) %>% 
  replace_na(list(U_mL = 0)) %>%
  filter(Genus_og == "Mastomys" | Genus_og == "Rattus", Protein =="NP") %>%
  add_column("Outlier" = FALSE) %>%
  mutate(Positive_NP = as.character(Positive_NP))


# Point edits for two rattus outlier points and filter out of sets
analysis2[1,47] = TRUE
analysis2[3,47] = TRUE

analysis2 <- analysis2 %>% mutate(Positive_GP = as.character(Positive_GP),
                     Positive_NP = as.character(Positive_NP))
analysis2[1,26] = "Outlier"
analysis2[3,27] = "Outlier"

analysisGP[3,47] = TRUE
analysisGP[3,26] = "Positive Outlier"
analysisGP <- analysisGP %>% 
              mutate(Positive_GP = fct_relevel(Positive_GP, c("FALSE", "TRUE", "Positive Outlier")))


analysisNP[9,47] = TRUE
analysisNP[9,27] = "Positive Outlier"
analysisNP <- analysisNP %>% 
              mutate(Positive_NP = fct_relevel(Positive_NP, c("FALSE", "TRUE", "Positive Outlier")))


analysisNP_outlier <- analysisNP %>% filter(Outlier == TRUE)
analysisNP_noout <- analysisNP %>% filter(Outlier == FALSE)

analysisGP_outlier <- analysisGP %>% filter(Outlier == TRUE)
analysisGP_noout <- analysisGP %>% filter(Outlier == FALSE)

analysis2_nooutliers <- analysis2 %>% filter(Outlier == FALSE)
analysis2_outliersline <- analysis2 %>% filter(Outlier == TRUE | Animal_ID == 4)

# Negative cutoffs (from experimental data)
Negatives <- data.frame(Genus_og = c("Mastomys", "Mastomys", "Rattus", "Rattus"), 
                        Protein = c("NP", "GP", "NP", "GP"),
                        Cutoff_U_mL = c(45.89, 28.20, 1.41, 5.24)) %>%   #from calculated cutoffs
                        mutate(Protein = fct_relevel(Protein, c("NP", "GP")))
NegativesGP <-  Negatives %>%
                filter(Protein == "GP")
NegativesNP <-  Negatives %>%
                filter(Protein == "NP")


# Shaprio-Wilk test for normality of distribution of differences for specimens with at least one positive IgG
analysis %>%
  filter(Genus_og == "Mastomys" | Genus_og == "Rattus") %>%
  filter(IgG_Positive == TRUE) %>%
  replace_na(list(U_mL_GP = 0, U_mL_NP = 0)) %>%
  mutate(difference = U_mL_NP - U_mL_GP) %>% 
  split(.$Genus_og) %>%
  map (~ shapiro_test(difference, data = .x))

### Calculating difference in mean IgG relative values for positive specimens

## Mastomys

# Mastomys paired specimens (specimens positive for both NP & GP)
paired_samples <- analysis2 %>%
  group_by(Genus_og, Animal_ID) %>%
  filter(n() == 2) %>%
  filter(Genus_og != "Rattus") %>%
  ungroup()

# Function to calculate mean difference between NP and GP (mastomys)
mean_diff <- function(data, indices) {
  d <- data[indices, ] 
  mean_NP <- mean(d$U_mL[d$Protein == "NP" & d$Positive_NP != "False"])
  mean_GP <- mean(d$U_mL[d$Protein == "GP" & d$Positive_GP != "False"])
  return(mean_NP - mean_GP)
}

# Mean U/mL for each protein, difference in means, and bootstrapped 95% CI of difference in means for positive mastomys specimen
paired_samples %>%
  filter(Genus_og == "Mastomys") %>%
  nest() %>%
  mutate(
    boot_res = map(data, ~boot(.x, mean_diff, R = 1000)),  
    CI = map(boot_res, ~boot.ci(.x, type = "perc")) 
  ) %>%
  mutate(
    lower_CI = map_dbl(CI, ~ifelse(is.null(.x$perc), NA, .x$perc[4])),
    upper_CI = map_dbl(CI, ~ifelse(is.null(.x$perc), NA, .x$perc[5])), 
    mean_difference = map_dbl(boot_res, ~mean(.x$t, na.rm = TRUE)),
    mean_GP = map_dbl(data, ~mean(.x$U_mL[.x$Positive_GP != "False" & .x$Protein == "GP"], na.rm = TRUE)),
    mean_NP = map_dbl(data, ~mean(.x$U_mL[.x$Positive_NP != "False" & .x$Protein == "NP"], na.rm = TRUE)) 
  ) %>%
  select(mean_NP, mean_GP, mean_difference, lower_CI, upper_CI)

# Paired Wilcox test
paired_tests <- paired_samples %>%
  group_by(Genus_og) %>%
  wilcox_test(U_mL ~ Protein, paired = TRUE) %>%
  add_significance() %>%
  mutate(test_type = "paired")

## Rattus 

# Positive rattus samples (running all unpaired since only 1/10 is paired NP+/GP+) 
unpaired_samples <- analysis2 %>%
  filter(!(Animal_ID %in% paired_samples$Animal_ID), 
         Genus_og != "Mastomys")

# Function to calculate mean difference between NP and GP (rattus)
mean_diff_rattus <- function(data, indices) {
  d <- data[indices, ]
  mean_NP <- mean(d$U_mL[d$Protein == "NP" & d$Positive_NP != "False"])
  mean_GP <- mean(d$U_mL[d$Protein == "GP" & d$Positive_GP != "False"])
  return(mean_GP - mean_NP) 
}

# Mean U/mL for each protein, difference in means, and bootstrapped 95% CI of difference in means for positive rattus specimens
unpaired_samples %>%
  filter(Genus_og == "Rattus") %>%
  nest() %>%
  mutate(
    boot_res = map(data, ~boot(.x, mean_diff_rattus, R = 1000)),
    CI = map(boot_res, ~boot.ci(.x, type = "perc"))
  ) %>%
  mutate(
    lower_CI = map_dbl(CI, ~ifelse(is.null(.x$perc), NA, .x$perc[4])),
    upper_CI = map_dbl(CI, ~ifelse(is.null(.x$perc), NA, .x$perc[5])), 
    mean_difference = map_dbl(boot_res, ~mean(.x$t, na.rm = TRUE)),
    mean_GP = map_dbl(data, ~mean(.x$U_mL[.x$Positive_GP != "False" & .x$Protein == "GP"], na.rm = TRUE)),  
    mean_NP = map_dbl(data, ~mean(.x$U_mL[.x$Positive_NP != "False" & .x$Protein == "NP"], na.rm = TRUE)) 
  ) %>%
  select(mean_GP, mean_NP, mean_difference, lower_CI, upper_CI)

# Mann-Whitney U test (unpaired wilcox test)
unpaired_tests <- unpaired_samples %>%
  group_by(Genus_og) %>%
  filter(n() > 1) %>% 
  wilcox_test(U_mL ~ Protein, paired = FALSE) %>%
  add_significance() %>%
  mutate(test_type = "unpaired")

## Combine paired and unpaired test results for figure
Fig2_wilcox <- bind_rows(paired_tests, unpaired_tests) %>%
  mutate(y.position = 4, 
         group1 = "GP", 
         group2 = "NP",
         p = case_when(p < 0.0001 ~ "<0.0001",
                       TRUE ~ formatC(p, format = "g", digits = 2)))

Fig2_wilcox_use <- Fig2_wilcox %>%
  filter(Genus_og == "Mastomys" & test_type == "paired" | 
         Genus_og == "Rattus" & test_type == "unpaired")

## Graph
# if change outlier.shape in geom_boxplot from NA to 1, will see that only outliers are the highest NP/GP rattus values. 
# this information was used to reflect the two outlier points in the data set (as seen above) 
ggplot() +
  stat_summary(data = analysis2, aes(Protein, y = U_mL), 
               fun.data = mean_se, fun.args = list(mult = 1), 
               geom = "crossbar")+
  geom_point(data = analysisGP, aes(Protein, U_mL, fill = Positive_GP, group = Animal_ID, shape = Positive_GP), 
            position = position_dodge(0.35)) +
  geom_point(data = analysisNP, aes(Protein, U_mL, fill = Positive_NP, group = Animal_ID, shape = Positive_NP), 
            position = position_dodge(0.35)) +
  geom_line(data = analysis2, aes(Protein, U_mL, group = Animal_ID), 
            size = 0.2, position = position_dodge(0.35))+
  geom_segment(aes(x = 0.5, y = Cutoff_U_mL, xend = 1.5, yend = Cutoff_U_mL, linetype = "Negative Cutoff"), 
               NegativesNP, color = "darkred", size = 0.75) +
  geom_segment(aes(x = 1.5, y = Cutoff_U_mL, xend = 2.5, yend = Cutoff_U_mL, linetype = "Negative Cutoff"), 
               NegativesGP, color = "darkred", size = 0.75) +
  scale_y_log10(name = "Mean IgG Relative Value \n (U/mL)" , 
                guide = "prism_offset",
                label = scales::label_number(accuracy = 1),
                limits = c(1,10000)) +
  facet_wrap(~Genus_og) +
  scale_fill_grey(name = "IgG Status", 
                  labels = c("Negative", "Positive", "Positive Outlier")) + 
  scale_shape_manual("IgG Status", 
                     labels = c("Negative", "Positive" ,"Positive Outlier"), 
                     values = c(21, 21, 10), 
                     guide = guide_legend(override.aes = list(shape = c(21, 21, 10)))) +
  scale_linetype_manual("       ", 
                        values=c("Negative Cutoff"=1), 
                        guide = guide_legend(override.aes = list(shape = NA)))+
  theme_pubr() +
  theme(legend.position = "right",
        legend.margin=unit(0, "cm"), 
        strip.text = element_text(size = 12, face = "italic"),
        axis.text.y = element_text(size = 10, vjust = 0.2), 
        axis.text.x = element_text(size = 10),
        axis.title.y =element_text(size = 12), 
        axis.title.x.bottom = element_text(size = 12)) +
  add_pvalue(Fig2_wilcox_use, label = "p = {p}", tip.length = 0.009, bracket.nudge.y = -0.009) #+
#ggsave("Figure2.png", width = 2300, height = 1500, units = "px", dpi = 320)

# 272, 445, 477 may be NP IgG positive when you look at the distribution with the cutoff line...

### Figure 3: Mean IgG relative value vs RDT Score for Mastomys, by protein ##########################################################

# Calculate Spearman's rho and p-value for IgG vs. RDT score for each protein (Mastomys only)
spearman <- analysis2 %>%
  filter(Genus_og == "Mastomys") %>%
  group_by(Protein) %>%
  summarise(spearman = list(broom::tidy(cor.test(U_mL, RDT_Score, method = "spearman")))) %>%
  unnest(cols = c(spearman)) %>%
  group_by(Protein) %>%
  summarize(label = paste0("rho == ", ifelse(Protein == "NP", round(estimate, 2), round(estimate, 3))), label2 =  paste0("~italic(p) == ", round(p.value, 2)))


# 95% CI of Spearman's rho for each protein
analysis2 %>%
  filter(Genus_og == "Mastomys") %>%
  group_by(Protein) %>%
  nest() %>%
  mutate(ci = map(data, ~tidy(ci.spear(0.05, .x$U_mL, .x$RDT_Score)))) %>% 
  unnest(ci) #[,3] and [,4] are LB and UB of CI respectively
        
# Graph
analysis2 %>%
  filter(Genus_og == "Mastomys") %>%
  ggplot(label = spearman$label)+
  geom_point(aes(x = RDT_Score, y = U_mL, fill = Ag_RDT), shape = 21) +
  geom_smooth(aes(x = RDT_Score, y = U_mL),
              color = "black",
              method = "lm",
              se = TRUE,
              fullrange = TRUE) +
  facet_wrap(~factor(Protein, levels = c("NP", "GP"))) + 
  geom_text(data = spearman, aes(x= 5, y = 10000, label = label), 
                   parse = TRUE, 
                   vjust = "top", 
                   hjust = "right", 
                   color = "black") +
  geom_text(data = spearman, aes(x= 5, y = 6000, label = label2), 
            parse = TRUE, 
            vjust = "top", 
            hjust = "right", 
            color = "black") +
  labs(y = "Mean IgG Relative Value \n (U/mL)", 
       x = "RDT Score") +
  scale_y_log10(oob = scales::squish_infinite, 
                limits = c(10, 10000), 
                breaks = scales::breaks_log(n=5, base = 10), 
                guide = "prism_offset_minor", 
                labels = scales::label_number()) +
  scale_fill_grey(name = "Antigen Status", labels = c("Ag-", "Ag+") ) +
  theme_pubr() +
  theme(legend.direction = "vertical", legend.position = "right",
        axis.text.y = element_text(size = 10, vjust = 0.2), 
        axis.text.x = element_text(size = 10, vjust = 0.2), 
        axis.title.y =element_text(size = 12), 
        axis.title.x.bottom = element_text(size = 12)) #+
  #ggsave("Figure3.png", width = 2300, height = 1500, units = "px", dpi = 320)


### Figure A1: Correlation between GP IgG and NP IgG relative values #################################################################

# Dataframe for Figure A1
analysis_abpos <- analysis %>%
  filter(IgG_Positive == "TRUE", Genus_og == "Mastomys") %>%
  mutate(Positive = case_when(Positive_GP == TRUE & Positive_NP == FALSE ~ "GP+/NP-",
                               Positive_GP == FALSE & Positive_NP == TRUE ~ "GP-/NP+", 
                               TRUE ~ "GP+/NP+"))

# Shaprio-Wilk test for normality of distribution of differences for specimens with at least one positive IgG
analysis_abpos %>%
  replace_na(list(U_mL_GP = 0, U_mL_NP = 0)) %>%
  mutate(difference = U_mL_NP - U_mL_GP) %>% 
  split(.$Genus_og) %>%
  map (~ shapiro_test(difference, data = .x))


# Spearman's rho and p-value of correlation between GP IgG value and NP IgG value
spearman_sup <- analysis_abpos %>%
  replace_na(list(U_mL_GP = 0, U_mL_NP = 0)) %>% # assign 0 for negative tests without a calculated U/mL
  summarise(spearman = list(broom::tidy(cor.test(U_mL_GP, U_mL_NP, method = "spearman", conf.level = 0.95, exact = TRUE)))) %>%
  unnest(cols = c(spearman)) %>%
  summarize(label = paste0("rho == ", round(estimate, 2)), 
            label2 = paste0("\U1D631 ", "<0.0001") )

# 95% CI of Spearman's rho 
analysis_abpos %>%
  filter(Genus_og == "Mastomys") %>%
  replace_na(list(U_mL_GP = 0, U_mL_NP = 0)) %>% 
  nest() %>%
  mutate(ci = map(data, ~tidy(ci.spear(0.05, .x$U_mL_GP, .x$U_mL_NP)))) %>% 
  unnest(ci)

# Graph  
analysis_abpos %>%
  ggplot()+
  geom_point(aes(U_mL_NP, U_mL_GP, fill = Positive), shape = 21) +
  geom_smooth(aes(U_mL_NP, U_mL_GP), 
              method = lm, 
              se = TRUE,
              fullrange = TRUE,
              color = "black",
              fill = "grey",
              alpha = 0.3)+
  geom_text(aes(x = Inf, y = Inf, label = label),
            parse = TRUE,
            vjust = "top", 
            hjust = "right", 
            color = "black", 
            data = spearman_sup) +
  geom_text(data = spearman_sup, aes(x= 9000, y = 5000, label = label2), 
            parse = FALSE, 
            vjust = "top", 
            hjust = "right", 
            color = "black") + 
  labs(y = "Mean GP IgG Relative Value \n (U/mL)", 
       x = "Mean NP IgG Relative Value \n (U/mL)") +
  scale_x_log10(oob = scales::squish_infinite, 
                limits = c(4e0, 9e3), 
                breaks = scales::breaks_log(n=5, base = 10), 
                guide = "prism_offset_minor", 
                labels = scales::label_number()) +
  scale_y_log10(oob = scales::squish_infinite, 
                limits = c(4e0, 9e3), 
                breaks = scales::breaks_log(n=5, base = 10), 
                guide = "prism_offset_minor", 
                labels = scales::label_number()) + 
  scale_fill_grey("IgG Status",
                  start = 1,
                  end = 0) +
  theme_pubr() +
  theme(legend.direction = "vertical", legend.position = "right",
        axis.text.y = element_text(size = 10, vjust = 0.2), 
        axis.text.x = element_text(size = 10, vjust = 0.2), 
        axis.title.y =element_text(size = 12), 
        axis.title.x.bottom = element_text(size = 12)) #+
#ggsave("FigureA1.png", width = 2300, height = 1500, units = "px", dpi = 320)
