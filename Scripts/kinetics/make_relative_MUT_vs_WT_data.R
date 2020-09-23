# this code takes the NMR< and GAP and GEF kinetics data and 
### a) calculates the relative value (relative to WT), including the propagated standard error of relative mean
### b) also adds the natural logarithm values with the appropriate standard error for ln-transformed (as se/mean)
library(tidyverse)

nmr_data <- read_csv('Data/31P_NMR_Data.csv', col_types = cols()) %>%
  filter(peak %in% c('gamma1','gamma2')) %>%
  select(variant, peak, integral_rel) %>%
  spread(peak, integral_rel) %>%
  mutate('gamma1' = case_when(is.na(gamma1) ~ 0.03*gamma2, T ~ gamma1),
         'gamma2' = case_when(is.na(gamma2) ~ 0.03*gamma1, T ~ gamma2)) %>%
  mutate('percent_gamma_NMR' = 100*gamma2/(gamma1+gamma2)) %>%
  mutate('gamma_ratio' = gamma2/gamma1) %>% 
  mutate('se' = NA) %>% 
  pivot_longer(c(percent_gamma_NMR, gamma_ratio), names_to = 'measure') %>% 
  select('mutant' = variant, se, value, measure) %>% 
  arrange(measure, mutant)

# Read in intrinsic and GAP-mediated GTP hydrolysis data
intrinsic_hydrolysis <- read_tsv('Data/RanGAP_assay/intrinsic_hydrolysis.txt') %>% 
  select(mutant, 'int' = mean_rel_rate, 'se' = se_rel_rate) %>% 
  gather("measure", "value", -mutant, -se) %>% 
  unique()

intrinsic_hydrolysis %>% 
  select(mutant, mean_rel_rate, se_rel_rate) %>% 
  unique() %>% write_tsv(., path = 'Revisions/Supplementary_Files/Supplementary_File_1/excel files/intrinsic_raw.txt')

GAP.effic <- read_tsv("Data/RanGAP_assay/GAP_kinetics_MichaelisMenten_parameters.txt") %>% 
  select(mutant, 'GAP_kcat_Km' = mean_kcat_Km, se) %>% 
  gather("measure", "value", -mutant, -se) %>% 
  unique()
GAP.kcat <- read_tsv("Data/RanGAP_assay/GAP_kinetics_MichaelisMenten_parameters.txt") %>% 
  select(mutant, 'GAP_kcat' = mean_kcat, 'se' = kcat_se) %>% 
  gather("measure", "value", -mutant, -se) %>% 
  unique()
GAP.Km <- read_tsv("Data/RanGAP_assay/GAP_kinetics_MichaelisMenten_parameters.txt") %>% 
  select(mutant, 'GAP_Km' = mean_Km, 'se' = Km_se) %>% 
  gather("measure", "value", -mutant, -se) %>% 
  unique()

GEF.effic <- read_tsv('Data/RanGEF_assay/GEF_kinetics_MichaelisMenten_parameters.txt') %>% 
  select(mutant, 'GEF_kcat_Km' = kcat_Km, 'se' = kcat_Km_se) %>% 
  gather("measure", "value", -mutant, -se) %>% 
  unique()
GEF.kcat <- read_tsv('Data/RanGEF_assay/GEF_kinetics_MichaelisMenten_parameters.txt') %>% 
  select(mutant, 'GEF_kcat' = kcat, 'se' = kcat_se) %>% 
  gather("measure", "value", -mutant, -se) %>% 
  unique()
GEF.Km <- read_tsv('Data/RanGEF_assay/GEF_kinetics_MichaelisMenten_parameters.txt') %>% 
  select(mutant, 'GEF_Km' = Km, 'se' = Km_se) %>% 
  gather("measure", "value", -mutant, -se) %>% 
  unique()
data <- bind_rows(GAP.effic, GAP.Km, GAP.kcat, intrinsic_hydrolysis, nmr_data, GEF.effic, GEF.Km, GEF.kcat)

WT_parameters <- data %>% 
  filter(mutant == "WT") %>% 
  mutate('value'= ifelse(measure == 'NMR', 1, value)) %>% 
  select(-mutant)

rel_data <- data %>% 
  inner_join(., WT_parameters, by = 'measure') %>% 
  select(mutant, measure, 'value' = value.x, 'se' = se.x, 'wt_value' = value.y, 'wt_se' = se.y) %>% 
  mutate('rel_value' = value/wt_value, 'rel_se' = (sqrt( (se/value)^2 + (wt_se/wt_value)^2 ) * value/wt_value)) %>% 
  select(mutant, measure, 'rel_to_WT' = rel_value, 'se' = rel_se) %>% 
  mutate('ln_rel_to_WT' = log(rel_to_WT), 'ln_se' = se/rel_to_WT)

### now get GAP/GEF values
# kcat/Km
gap_temp <- rel_data %>% 
  filter(measure == 'GAP_kcat_Km') %>% 
  select(mutant, 'GAP_rel_to_WT' = rel_to_WT, 'GAP_se' = se)
temp <- rel_data %>% 
  filter(measure == 'GEF_kcat_Km') %>% 
  select(mutant, 'GEF_rel_to_WT' = rel_to_WT, 'GEF_se' = se) %>% 
  inner_join(., gap_temp, by = 'mutant') %>% 
  mutate('GAP/GEF kcat/Km' = GAP_rel_to_WT/GEF_rel_to_WT) %>% 
  mutate('GAP/GEF kcat/Km se' = (sqrt( (GAP_se/GAP_rel_to_WT)^2 + (GEF_se/GEF_rel_to_WT)^2 ) * GAP_rel_to_WT/GEF_rel_to_WT)) %>% 
  select(mutant, `GAP/GEF kcat/Km`, 'GAP/GEF kcat/Km se') %>% 
  mutate('measure' = 'GAP/GEF kcat/Km', 'rel_to_WT' = `GAP/GEF kcat/Km`, se = `GAP/GEF kcat/Km se`) %>% 
  # mutate('rel_to_WT' = ifelse(rel_to_WT > 40, rel_to_WT/7, rel_to_WT)) %>% 
  # mutate('se' = ifelse(rel_to_WT > 40, se/7, se)) %>% 
  mutate('ln_rel_to_WT' = log(rel_to_WT), 'ln_se' = se/rel_to_WT) %>% 
  select(mutant, measure, rel_to_WT, se, ln_rel_to_WT, ln_se)
rel_data <- bind_rows(rel_data, temp)
  
# kcat
gap_temp <- rel_data %>% 
  filter(measure == 'GAP_kcat') %>% 
  select(mutant, 'GAP_rel_to_WT' = rel_to_WT, 'GAP_se' = se)
temp <- rel_data %>% 
  filter(measure == 'GEF_kcat') %>% 
  select(mutant, 'GEF_rel_to_WT' = rel_to_WT, 'GEF_se' = se) %>% 
  inner_join(., gap_temp, by = 'mutant') %>% 
  mutate('GAP/GEF kcat' = GAP_rel_to_WT/GEF_rel_to_WT) %>% 
  mutate('GAP/GEF kcat se' = (sqrt( (GAP_se/GAP_rel_to_WT)^2 + (GEF_se/GEF_rel_to_WT)^2 ) * GAP_rel_to_WT/GEF_rel_to_WT)) %>% 
  select(mutant, `GAP/GEF kcat`, 'GAP/GEF kcat se') %>% 
  mutate('measure' = 'GAP/GEF kcat', 'rel_to_WT' = `GAP/GEF kcat`, se = `GAP/GEF kcat se`) %>% 
  mutate('ln_rel_to_WT' = log(rel_to_WT), 'ln_se' = se/rel_to_WT) %>% 
  select(mutant, measure, rel_to_WT, se, ln_rel_to_WT, ln_se)
rel_data <- bind_rows(rel_data, temp)
  
# Km
gap_temp <- rel_data %>% 
  filter(measure == 'GAP_Km') %>% 
  select(mutant, 'GAP_rel_to_WT' = rel_to_WT, 'GAP_se' = se)
temp <- rel_data %>% 
  filter(measure == 'GEF_Km') %>% 
  select(mutant, 'GEF_rel_to_WT' = rel_to_WT, 'GEF_se' = se) %>% 
  inner_join(., gap_temp, by = 'mutant') %>% 
  mutate('GAP/GEF Km' = GAP_rel_to_WT/GEF_rel_to_WT) %>% 
  mutate('GAP/GEF Km se' = (sqrt( (GAP_se/GAP_rel_to_WT)^2 + (GEF_se/GEF_rel_to_WT)^2 ) * GAP_rel_to_WT/GEF_rel_to_WT)) %>% 
  select(mutant, `GAP/GEF Km`, 'GAP/GEF Km se') %>% 
  mutate('measure' = 'GAP/GEF Km', 'rel_to_WT' = `GAP/GEF Km`, se = `GAP/GEF Km se`) %>% 
  mutate('ln_rel_to_WT' = log(rel_to_WT), 'ln_se' = se/rel_to_WT) %>% 
  select(mutant, measure, rel_to_WT, se, ln_rel_to_WT, ln_se)
rel_data <- bind_rows(rel_data, temp)

write_tsv(rel_data, path = 'Data/kinetics_data_relative_to_WT.txt')
