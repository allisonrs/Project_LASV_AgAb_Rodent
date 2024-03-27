library(tidyverse)
library(sf)
library(janitor)
library(gtsummary)
library(gt)
library(rstatix)
library(ggpubr)
library(ggprism)
library(tayloRswift)
library(epiR)
library(statpsych)
library(exact2x2)

### Prepare Data ##############################################################################

original_data <- read_csv("data/LASV_Analysis.csv") 
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

### Tables and Figures in code ###############################

# Table 1: Counts and row percentage of Ag+
# Table 2: Counts and row percentage of NP IgG+
# Table 3: Counts and row percentage of GP IgG+, NP IgG+, GP and NP IgG+, IgG+
# Table 4: Contingency table of Ag/IgG +/- for Mastomys
# Table S1: GPS Coordinates of All Villages in Study
# Table S2: Specimen counts by Cytochrome b PCR
# Table S3: Contingency Table of Ag/qRT-PCR status of serum of Ag+/- specimens
# Table S4: Sensitivity and specificity of serum Ag to qRT-PCR. 
# Table S5: qRT-PCR status of select organs from select Ag+ specimens 
# Table S6: Counts by Genus ID in field; mouse/rat reagents used
# Table S7: Contingency table of Ag/IgG +/- for Rattus
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
# serum counts
analysis %>%
  filter(!is.na(qPCR_Positive)) %>%
  drop_na(Ag_RDT, Positive_AVL) %>%
  count(Positive_AVL)

# Serum specificty/sensitivity
freqtable_RDTPCR_serum <- analysis %>%
  filter(!is.na(qPCR_Positive)) %>%
  drop_na(Ag_RDT, Positive_AVL) %>% 
  mutate(Ag_RDT = if_else(Ag_RDT == TRUE, "Ag+", "Ag-"),
         qRTPCR = if_else(Positive_AVL == TRUE, "qRT-PCR+", "qRT-PCR-")) %>%
  mutate(Ag_RDT = factor(Ag_RDT, levels = c("Ag+", "Ag-")),
         qRTPCR = factor(qRTPCR, levels = c("qRT-PCR+", "qRT-PCR-"))) %>%
  xtabs(~Ag_RDT+qRTPCR, data = .) 

## mcnemar_test
freqtable_RDTPCR_serum %>% mcnemar.exact()

### Table S4 ############################################################
serum_SS_RDTPCR_gt <- epi.tests(freqtable_RDTPCR_serum) %>%
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

### Table S5 ############################################################
## organs only, no serum
PCR_tissues <- analysis %>%
  filter(!is.na(qPCR_Positive), Ag_RDT == T) %>% 
  mutate(AVL_only = if_else(!is.na(Positive_AVL) 
                            & is.na(Positive_Liver)
                            & is.na(Positive_Lung) 
                            & is.na(Positive_Kidney) 
                            &is.na(Positive_Spleen), 
                            "TRUE", "FALSE")) %>% 
  filter(AVL_only == F)

PCR_tissues %>%
  group_by(qPCR_Positive) %>%
  count()

## organ counts
PCR_tissues_counts <- PCR_tissues %>%
  rename_with(~gsub("Positive_", "", .x)) %>%
  pivot_longer(cols = c(AVL, Liver, Lung, Kidney, Spleen), names_to ="organ", values_to = "PCR_Positive_organ", values_drop_na = TRUE) %>% 
  filter(organ != "AVL") 

PCR_tissues_counts%>%
  group_by(organ) %>%
  count(PCR_Positive_organ) %>%
  pivot_wider(names_from = PCR_Positive_organ, values_from = n)


### Table S6 ############################################################
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

# dataframes for Figure 2
analysis2 <- analysis %>% ## only IgG positive for at least one protein
  filter(IgG_Positive == "TRUE") %>%
  pivot_longer(cols = c("U_mL_GP", "U_mL_NP"), names_to = "Protein", values_to = "U_mL") %>% 
  mutate(Protein = if_else(Protein == "U_mL_GP", "GP", "NP")) %>%
  filter(Positive_GP == "TRUE" & Protein == "GP" | Positive_NP == "TRUE" & Protein == "NP") %>%
  mutate(Protein = fct_relevel(Protein, c("NP", "GP"))) %>%
  add_column(Outlier = FALSE)

analysisGP <- analysis %>%
  pivot_longer(cols = c("U_mL_GP", "U_mL_NP"), names_to = "Protein", values_to = "U_mL") %>% 
  mutate(Protein = if_else(Protein == "U_mL_GP", "GP", "NP")) %>% 
  replace_na(list(U_mL = 0)) %>%
  filter(Genus_og == "Mastomys" | Genus_og == "Rattus", Protein =="GP") %>%
  add_column("Outlier" = FALSE) %>%
  analysisNP %>% mutate(Positive_GP = as.character(Positive_GP))
  

analysisNP <- analysis %>%
  pivot_longer(cols = c("U_mL_GP", "U_mL_NP"), names_to = "Protein", values_to = "U_mL") %>% 
  mutate(Protein = if_else(Protein == "U_mL_GP", "GP", "NP")) %>% 
  replace_na(list(U_mL = 0)) %>%
  filter(Genus_og == "Mastomys" | Genus_og == "Rattus", Protein =="NP") %>%
  add_column("Outlier" = FALSE) %>%
  analysisNP %>% mutate(Positive_NP = as.character(Positive_NP))


# point edits for two rattus outlier points and filter out of sets
analysis2[1,47] = TRUE
analysis2[3,47] = TRUE

analysis2 <- analysis2 %>% mutate(Positive_GP = as.character(Positive_GP),
                     Positive_NP = as.character(Positive_NP))
analysis2[1,26] = "Outlier"
analysis2[3,27] = "Outlier"

analysisGP[3,47] = TRUE
analysisGP[3,26] = "Positive Outlier"
analysisGP <- analysisGP %>% mutate(Positive_GP = fct_relevel(Positive_GP, c("FALSE", "TRUE", "Positive Outlier")))


analysisNP[9,47] = TRUE
analysisNP[9,27] = "Positive Outlier"
analysisNP <- analysisNP %>% mutate(Positive_NP = fct_relevel(Positive_NP, c("FALSE", "TRUE", "Positive Outlier")))



analysiNP_outlier <- analysisNP %>% filter(Outlier == TRUE)
analysisNP_noout <- analysisNP %>% filter(Outlier == FALSE)

analysisGP_outlier <- analysisGP %>% filter(Outlier == TRUE)
analysisGP_noout <- analysisGP %>% filter(Outlier == FALSE)

analysis2_nooutliers <- analysis2 %>% filter(Outlier == FALSE)
analysis2_outliersline <- analysis2 %>% filter(Outlier == TRUE | Animal_ID == 4)

# Negatives 
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

# means of positives
analysis2 %>%
  group_by(Genus_og, Protein) %>%
  summarize(mean = mean(U_mL)) %>%
  as.data.frame()

# median of positives
analysis2 %>%
  group_by(Genus_og, Protein) %>%
  summarize(median = median(U_mL)) %>%
  as.data.frame()

# two-tailed Wilcox signed rank test [paired] of mean GP-Specific IgG vs. mean NP-specific IgG for Positives 
Fig2_wilcox <- analysis2 %>% 
  group_by(Genus_og) %>%
  nest() %>% 
  mutate(wilcox = map(data, ~wilcox_test(U_mL ~ Protein, data = .x), paired= TRUE)) %>%
  unnest(wilcox) %>%
  add_significance() %>%
  mutate(y.position = 4, 
         group1 = "GP", 
         group2 ="NP") %>%
  mutate(p = if_else(p.signif ==  "ns", 0.16, 0.00016))

### Graph 
# if change outlier.shape in geom_boxplot from NA to 1, will see that only outliers are the highest NP/GP rattus values. 
# this information was used to reflect the two outlier points in the data set (as seen above) 
analysis2 %>% 
  ggplot() +
  geom_boxplot(aes(Protein, (U_mL)), 
               varwidth = TRUE, outlier.shape = NA) + 
  geom_point(aes(Protein, (U_mL), fill = Positive_GP, group = Animal_ID, shape = Positive_GP), 
             analysisGP, position = position_dodge(0.35)) +
  geom_point(aes(Protein, (U_mL), fill = Positive_NP, group = Animal_ID, shape = Positive_NP), 
             analysisNP, position = position_dodge(0.35)) +
  geom_line(aes(Protein, (U_mL), group = Animal_ID), 
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
  scale_x_discrete(name = "Protein") +
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
  add_pvalue(Fig2_wilcox, label = "p = {p}", tip.length = 0.009, bracket.nudge.y = -0.009) #+
  #ggsave("Figure2_edit5.png", width = 2300, height = 1500, units = "px", dpi = 320)

# 272, 445, 477 may be positive when you look at the distribution with the cutoff line...

### Figure 3 ####################################################################

#spearman correlation 
spearman <- analysis2 %>%
  filter(Genus_og == "Mastomys") %>%
  group_by(Protein) %>%
  summarise(spearman = list(broom::tidy(cor.test(U_mL, RDT_Score, method = "spearman")))) %>%
  unnest(cols = c(spearman)) %>%
  group_by(Protein) %>%
  summarize(label = paste0("rho == ", ifelse(Protein == "NP", round(estimate, 2), round(estimate, 3))), label2 =  paste0("~italic(p) == ", round(p.value, 2)))


# 95% CI for spearman's correlation 
analysis2 %>%
  filter(Genus_og == "Mastomys") %>%
  group_by(Protein) %>%
  nest() %>%
  mutate(ci = map(data, ~tidy(ci.spear(0.05, .x$U_mL, .x$RDT_Score)))) %>% 
  unnest(ci) ## estimate, SE, LL, UL
        
#double check results of mapped
analysis2GP <- analysis2 %>% filter(Genus_og == "Mastomys", Protein == "GP") 
ci.spear(0.05, analysis2GP$U_mL, analysis2GP$RDT_Score)

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
  scale_fill_taylor(name = "Antigen Status", palette = "reputation", labels = c("Negative", "Positive") ) +
  theme_pubr() +
  theme(legend.direction = "vertical", legend.position = "right",
        axis.text.y = element_text(size = 10, vjust = 0.2), 
        axis.text.x = element_text(size = 10, vjust = 0.2), 
        axis.title.y =element_text(size = 12), 
        axis.title.x.bottom = element_text(size = 12)) #+
#ggsave("Figure3_Edits4.png", width = 2300, height = 1500, units = "px", dpi = 320)


### Table 3 (Mastomys Contingency Table) and Table S7 (Rattus Contingency Table) ###############
  
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
  addmargins() 

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




### Figure S1 ###########################################################################

# Dataframe for figure S1
analysis_abpos <- analysis %>%
  filter(IgG_Positive == "TRUE", Genus_og == "Mastomys") 

# Shaprio-Wilk test for normality of distribution of differences for specimens with at least one positive IgG
analysis_abpos %>%
  replace_na(list(U_mL_GP = 0, U_mL_NP = 0)) %>%
  mutate(difference = U_mL_NP - U_mL_GP) %>% 
  split(.$Genus_og) %>%
  map (~ shapiro_test(difference, data = .x))



#spearman's rho and p value for figure S1
spearman_sup <- analysis_abpos %>%
  summarise(spearman = list(broom::tidy(cor.test(U_mL_GP, U_mL_NP, method = "spearman")))) %>%
  unnest(cols = c(spearman)) %>%
  summarize(label = paste0("rho == ", round(estimate, 2)), 
            label2 = paste0("\U1D631", " < 0.0001"))
#spearman 95% CI
ci.spear(0.05, analysis_abpos$U_mL_GP, analysis_abpos$U_mL_NP)

# Graph  
analysis_abpos %>%
  ggplot()+
  geom_point(aes(U_mL_NP, U_mL_GP)) +
  geom_smooth(aes(U_mL_NP, U_mL_GP), 
              method = lm, 
              se = TRUE,
              fullrange = TRUE,
              color = "black")+
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
  scale_fill_taylor(palette = "reputation")+
  theme_pubr() +
  theme(legend.direction = "vertical", legend.position = "right",
        axis.text.y = element_text(size = 10, vjust = 0.2), 
        axis.text.x = element_text(size = 10, vjust = 0.2), 
        axis.title.y =element_text(size = 12), 
        axis.title.x.bottom = element_text(size = 12)) #+
#ggsave("FigureS1_edits3.png", width = 2300, height = 1500, units = "px", dpi = 320)
