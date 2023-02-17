library(tidyverse)
library(sf)
library(janitor)
library(gtsummary)
library(rstatix)
library(ggpubr)
library(ggprism)
library(tayloRswift)

### Prepare Data ##############################################################################

original_data <- read_csv("data/LASV_Analysis.csv") ## output of Interpolated_Ab_Concentration.R file
analysis <- original_data %>%
                    filter(!is.na(Ag_RDT), !is.na(IgG_Positive)) ## excludes specimens not tested



### Shape Files -- data source for SL Maps: https://data.humdata.org/dataset/cod-ab-sle ##########

# Shape File of SL Provinces
sle_adm1 <- read_sf("data/SL_Map_Data/sle_admbnda_adm1_1m_gov_ocha_20161017.shp") %>%   
  clean_names()

# Shape File of SL Districts
sle_adm2 <- read_sf("data/SL_Map_Data/sle_admbnda_adm2_1m_gov_ocha_20161017.shp") %>% 
  clean_names() %>%
  filter(admin2name %in% "Kenema") ## filter to include Kenema district only

### Tables and Figures in code are listed in order of generation ###############################

# Table 1: Counts and row percentage of Ag+
# Table 2: Counts and row percentage of NP IgG+
# Table 3: Counts and row percentage of GP IgG+, NP IgG+, GP and NP IgG+, IgG+
# Table 4: Contingency table of Ag/IgG +/- for Mastomys
# Table 5: Contingency table of Ag/IgG +/- for Rattus
# Table S1: GPS Coordinates of All Villages in Study
# Table S2: Specimen counts by Cytochrome b PCR
# Table S3: Counts by Genus ID in field; mouse/rat reagents used
# Figure 1: Location of villages sampled within provincial level map of Sierra Leone
# Figure 2: Distribution of GP and NP IgG by protein and genus
# Figure 3: Mean IgG concentration vs antigen strength for Mastomys
# Figure S1: Correlation between strength of GP IgG vs NP IgG in Mastomys


## Table S1 ############################################################
GPS_Coordinates_Villages <- original_data %>%
  group_by(Village) %>%
  summarize(first(Lat_Convert), first(Long_Use)) %>%
  rename(lat = `first(Lat_Convert)`, lon = `first(Long_Use)`) 


### Table S2 ############################################################
analysis %>%
  group_by(Cytb_species) %>%
  count() %>%
  arrange(desc(n)) %>%
  replace_na(list(Cytb_species = "N/A: Not tested")) 


### Table S3 ############################################################
analysis %>%
  group_by(Genus_og) %>%
  count() %>%
  arrange(desc(n)) %>%
  rename(Field_ID_Genus = Genus_og) %>%
  mutate(Ab_IgG_ELISA = if_else(Field_ID_Genus == "Rattus", "Rat", "Mouse"), .before = n) 


### Figure 1, using GPS coordinates generated in Table S1 ###############
ggplot() + 
  geom_sf(data = sle_adm1, fill = NA) + 
  geom_sf(data = sle_adm2, mapping = aes(fill = admin2name), fill = "gray") + 
  geom_point(data = GPS_Coordinates_Villages, mapping = aes(x = lon, y = lat), shape = 1, size = 0.3) +
  coord_sf() +
  labs(x = "Longitude", y = "Latitude") +
  theme(panel.background = element_rect(fill = 'white', colour = 'white'), 
        legend.key = element_rect(fill = 'white', color = 'white'), 
        legend.background = element_rect(fill = "white", color = "white")) #+
  #ggsave("Figure1.png", width = 1500, height = 1500, units = "px", dpi = 320) 


### correct data frame for the one Malacomys later determined to be a Praomys ##########
analysis[6,10] <- "Praomys" 


### Table 1 ############################################################################
analysis %>%
  select(Genus_og, Ag_RDT) %>%
  rename("Genus" = "Genus_og") %>%
  tbl_summary(by = Ag_RDT, sort = list(everything() ~ "frequency"), type = list(everything() ~ "categorical"), percent = "row", missing = "no") %>%
  add_overall() 


### Table 2 ############################################################################
analysis %>%
  select(Genus_og, Positive_NP) %>%
  rename("Genus" = "Genus_og") %>%
  tbl_summary(by = Positive_NP, sort = list(everything() ~ "frequency"), type = list(everything() ~ "categorical"), percent = "row", missing = "no") %>%
  add_overall() 


### Table 3 ############################################################################
analysis %>%
  mutate(Genus_og = factor(Genus_og, levels = c("Mastomys", "Rattus", "Praomys", "Hylomyscus"))) %>%
  mutate(IgG = case_when(Positive_GP == TRUE & Positive_NP == FALSE ~ "GP IgG+ only",
                         Positive_GP == FALSE & Positive_NP == TRUE ~ "NP IgG+ only",
                         Positive_GP == TRUE & Positive_NP == TRUE ~ "GP and NP IgG+",
                         TRUE ~ "IgG-")) %>%
  select(Genus_og, IgG, IgG_Positive) %>%
  rename("Genus" = "Genus_og") %>%
  tbl_summary(by = Genus, sort = list(everything() ~ "frequency"), type = list(everything() ~ "categorical"), percent = "column", missing = "no") %>%
  add_overall(last = TRUE) %>%
  add_n()

### Figure 2 ###########################################################################

# dataframe for Figure 2
analysis2 <- analysis %>%
  filter(IgG_Positive == "TRUE") %>%
  pivot_longer(cols = c("U_mL_GP", "U_mL_NP"), names_to = "Protein", values_to = "U_mL") %>% 
  mutate(Protein = if_else(Protein == "U_mL_GP", "GP", "NP")) %>%
  filter(Positive_GP == "TRUE" & Protein == "GP" | Positive_NP == "TRUE" & Protein == "NP")
  

# unpaired two-tailed Student's t-test
Fig2_ttest <- analysis2 %>% 
  group_by(Genus_og) %>%
  nest() %>% 
  mutate(t_test = map(data, ~t_test(U_mL ~ Protein, data = .x))) %>%
  unnest(t_test) %>%
  add_significance() %>%
  mutate(y.position = ifelse(p.signif == "ns", NA, 4))

# Graph with t-test results 
analysis2 %>%
  ggplot() +
  geom_boxplot(aes(Protein, log10(U_mL)), varwidth = TRUE, outlier.shape = NA) +
  geom_point(aes(Protein, log10(U_mL), fill = Protein, group = Animal_ID), position = position_dodge(0.35), shape = 21) +
  geom_line(aes(Protein, log10(U_mL), group = Animal_ID), size = 0.2, position = position_dodge(0.35))+
  coord_cartesian(ylim = c(0,4))+
  scale_y_continuous(name = "Mean IgG Concentration \n (U/mL)" , 
                    guide = "prism_offset",
                    label = c("1+e00", "1+e01", "1+e02", "1+e03", "1+e04")) +
  facet_wrap(~Genus_og, )+
  scale_x_discrete(name = "Protein") +
  scale_fill_taylor(palette = "reputation")+
  theme_pubr() +
  theme(legend.position = "right",
        strip.text = element_text(size = 12, face = "italic"),
        axis.text.y = element_text(size = 10, vjust = 0.2), 
        axis.text.x = element_text(size = 10),
        axis.title.y =element_text(size = 12), 
        axis.title.x.bottom = element_text(size = 12)) +
  add_pvalue(Fig2_ttest, label = "p.signif") #+
  #ggsave("Figure2.png", width = 2300, height = 1500, units = "px", dpi = 320)

### Figure 3 ####################################################################

# Linear regression of Ab concentration vs. RDT Score by Protein in Mastomys
Fig3_linreg <- analysis2 %>%
  filter(Ag_RDT == "TRUE", Genus_og == "Mastomys") %>%
  group_by(Protein) %>%
  select(Protein, Animal_ID, U_mL, RDT_Score) %>%
  nest() %>%
  mutate(linreg = map(data, ~lm(.$U_mL ~ .$RDT_Score)), 
         tidied = map(linreg, tidy), 
         glanced = map(linreg, broom::glance)) %>%
  unnest(glanced) %>%
  select(Protein, r.squared)

# Graph with linear regression 
analysis2 %>%
  filter(Ag_RDT == "TRUE", Genus_og == "Mastomys") %>%
  group_by(Protein) %>%
  ggplot()+
  geom_point(aes(x = RDT_Score, y = U_mL, fill = Protein), shape = 21) +
  facet_wrap(~Protein, ) +
  stat_smooth(aes(x = RDT_Score, y = U_mL), 
              color = "black", 
              se = TRUE,
              method = lm, 
              fullrange = TRUE) +
  geom_text(aes(x = Inf, y = Inf, label = paste0("~R^{2} == ", round(r.squared, 4))), 
            parse = TRUE,
            vjust = "top", 
            hjust = "right", 
            color = "black", 
            data = Fig3_linreg) +
  labs(y = "Mean IgG Concentration \n (U/mL)", 
       x = "RDT Score") +
  scale_y_log10(oob = scales::squish_infinite, 
                limits = c(1e1, 1e4), 
                breaks = scales::breaks_log(n=5, base = 10), 
                guide = "prism_offset_minor", 
                labels = scales::label_scientific()) +
  scale_fill_taylor(palette = "reputation") +
  theme_pubr() +
  theme(legend.direction = "vertical", legend.position = "right",
        axis.text.y = element_text(size = 10, vjust = 0.2), 
        axis.text.x = element_text(size = 10, vjust = 0.2), 
        axis.title.y =element_text(size = 12), 
        axis.title.x.bottom = element_text(size = 12)) #+
  #ggsave("Figure3.png", width = 2300, height = 1500, units = "px", dpi = 320)

### Table 4 (Mastomys Contingency Table) and Table 5 (Rattus Contingency Table) ###############
  
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
  xtabs(~Ag_RDT+IgG_Positive, data = .) %>%
  addmargins()

# Contingency tables and p-values
analysis %>%
  filter(Genus_og == "Mastomys" | Genus_og == "Rattus") %>%
  mutate(Ag_RDT = if_else(Ag_RDT == TRUE, "Ag +", "Ag -"),
         IgG = if_else(IgG_Positive == TRUE, "IgG +", "IgG -")) %>%
  group_by(Genus_og) %>%
  select(Genus_og, Ag_RDT, IgG) %>%
  nest() %>%
  mutate(contingency = map(data, ~xtabs(~Ag_RDT+IgG, data = .)),
         fisher = map(contingency, fisher.test), 
         contingency = map(contingency, addmargins),
         glanced = map(fisher, broom::glance)) %>%
  #mutate(contingency = map(contingency, as.data.frame.matrix)) %>%
  unnest(c(glanced)) %>%
  select(Genus_og, contingency, p.value)

### Figure S1 ###########################################################################

# Dataframe for figure S1
analysis_abpos <- analysis %>%
  filter(IgG_Positive == "TRUE", Genus_og == "Mastomys")

# R squared value for figure S1
rsq_label_abpos <- broom::glance(lm(U_mL_GP ~ U_mL_NP, data = analysis_abpos))

# Graph with linear regression
analysis_abpos %>%
  ggplot() +
  geom_point(aes(U_mL_NP, U_mL_GP)) +
  stat_smooth(aes(U_mL_NP, U_mL_GP), 
              se = TRUE,
              method = lm, 
              fullrange = TRUE,
              color = "black")+
  geom_text(aes(x = Inf, y = Inf, label = paste0("~R^{2} == ", round(r.squared, 4))),
            parse = TRUE,
            vjust = "top", 
            hjust = "right", 
            color = "black", 
            data = rsq_label_abpos) +
  labs(y = "Mean GP IgG Concentration \n (U/mL)", 
       x = "Mean NP IgG Concentration \n (U/mL)") +
  scale_x_log10(oob = scales::squish_infinite, 
                limits = c(4e0, 9e3), 
                breaks = scales::breaks_log(n=5, base = 10), 
                guide = "prism_offset_minor", 
                labels = scales::label_scientific()) +
  scale_y_log10(oob = scales::squish_infinite, 
                limits = c(4e0, 9e3), 
                breaks = scales::breaks_log(n=5, base = 10), 
                guide = "prism_offset_minor", 
                labels = scales::label_scientific()) + 
  scale_fill_taylor(palette = "reputation")+
  theme_pubr() +
  theme(legend.direction = "vertical", legend.position = "right",
        axis.text.y = element_text(size = 10, vjust = 0.2), 
        axis.text.x = element_text(size = 10, vjust = 0.2), 
        axis.title.y =element_text(size = 12), 
        axis.title.x.bottom = element_text(size = 12)) #+
  #ggsave("FigureS1.png", width = 2300, height = 1500, units = "px", dpi = 320)
  