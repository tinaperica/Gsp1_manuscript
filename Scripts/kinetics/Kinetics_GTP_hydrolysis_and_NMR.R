library(tidyverse)
library(ggrepel)
source('ucsf_colors.R')
main_directory <- 'Figure3_Biophysics/Plots'
extended_directory <- 'Extended_Figures/'
supplemental_directory <- 'Supplemental_Figures/'

# Read in NMR data
spectra_to_include <- c('H141R','H141V','K132H','Y157A','D79S','WT',
                        'G80A','R108L','R78K','T34G','Q147E',
                        'T34D','T34L','T34A','T34Q','T34E')

nmr_data <- read_csv('Data/31P_NMR_Data.csv', col_types = cols()) %>%
  filter(peak %in% c('gamma1','gamma2')) %>%
  select(variant, peak, integral_rel) %>%
  spread(peak, integral_rel) %>%
  mutate('gamma1' = case_when(is.na(gamma1) ~ 0.0, T ~ gamma1),
         'gamma2' = case_when(is.na(gamma2) ~ 0.0, T ~ gamma2)) %>%
  mutate('value' = 100*gamma2/(gamma1+gamma2), 'measure' = 'NMR', 'sd' = NA) %>%
  select('mutant' = variant, sd, value, measure) %>% 
  filter(mutant %in% spectra_to_include)

# Read in intrinsic and GAP-mediated GTP hydrolysis data
intrinsic_hydrolysis <- read_tsv('Data/RanGAP_assay/intrinsic_hydrolysis.txt') %>% 
  select(mutant, 'int'= mean_log_rel_rate, 'sd' = sd_log_rel_rate) %>% 
  gather("measure", "value", -mutant, -sd)

MM.data <- read_tsv("Data/RanGAP_assay/GAP_kinetics_MichaelisMenten_parameters.txt") %>% 
  select(mutant, 'kcat_Km' = log_kcat_over_Km, 'sd' = log_sd) %>% 
  gather("measure", "value", -mutant, -sd)
data <- bind_rows(MM.data, intrinsic_hydrolysis, nmr_data)

WT_parameters <- data %>% 
  filter(mutant == "WT") %>% 
  mutate('value'= ifelse(measure == 'NMR', 1, value)) %>% 
  select(-mutant)

data <- data %>% 
  inner_join(., WT_parameters, by = 'measure') %>% 
  select(mutant, measure, 'value' = value.x, 'sd' = sd.x, 'wt_value' = value.y, 'wt_sd' = sd.y) %>% 
  #mutate('rel_value' = value/wt_value, 'rel_sd' = (sqrt( (sd/value)^2 + (wt_sd/wt_value)^2 ) * value/wt_value)) %>% 
  mutate('rel_value' = ifelse(measure == 'NMR',value/wt_value, value - wt_value), 'rel_sd' = sqrt( sd^2 + wt_sd^2)) %>% 
  select(mutant, measure, 'value' = rel_value, 'sd' = rel_sd)

# plot GAP vs state 2
GAP_NMR_spread <- data %>%
  select(-sd) %>%
  filter(measure %in% c('kcat_Km', 'NMR')) %>%
  spread(measure, value)
GAP_NMR_spread <- data %>% 
  select(mutant, measure, sd) %>% 
  filter(measure == 'kcat_Km') %>% 
  inner_join(., GAP_NMR_spread, by = 'mutant') %>% 
  filter(! is.na(NMR))
outliers <- GAP_NMR_spread %>% 
  filter(mutant %in% c('R78K', 'D79S', 'K132H')) %>% 
  #mutate('kcat_Km' = log(kcat_Km)) %>% 
  summarise('ymin' = min(kcat_Km) - abs(0.1 * min(kcat_Km)),
            'ymax' = max(kcat_Km) + abs(0.1 * max(kcat_Km)),
            'xmin' = min(NMR) - abs(0.1 * min(NMR)),
            'xmax' = max(NMR) + abs(0.1 * max(NMR)))
filtered_GAP_NMR_spread <- GAP_NMR_spread %>% 
  filter(! mutant %in% c('R78K', 'D79S', 'K132H'))
GAPNMRfit <- lm(kcat_Km ~ NMR, data = filtered_GAP_NMR_spread)
GAP_NMR_spread %>% 
  ggplot(aes(y = kcat_Km, x = NMR)) + 
  geom_point(size = 2, alpha = 0.5) +
  geom_errorbar(aes(ymin = kcat_Km - sd, ymax = kcat_Km + sd), width = 2, size = 0.1) + 
  geom_rect(data = outliers, inherit.aes = F, fill = ucsf_colors$gray3, alpha = 0.2,
            aes(xmin = outliers$xmin, xmax = outliers$xmax, ymin = outliers$ymin, ymax = outliers$ymax)) + 
  geom_text_repel(aes(label = mutant), size = 2, point.padding = 0.1, segment.size = 0.3) +
  ylab(expression("GAP mediated GTP hydrolysis relative k"['cat']*"/K"['m'])) +
  xlab(expression("\n% "*gamma*" phosphate in "*gamma*" state 2")) +
  theme_classic() +
  geom_abline(intercept = GAPNMRfit$coefficients[1], slope = GAPNMRfit$coefficients[2], color = ucsf_colors$pink1, alpha = 0.5) +  
  theme(text = element_text(size = 6),
        axis.title.x = element_text(size = 6), 
        axis.title.y = element_text(size = 6),
        axis.line = element_line(size = 0.1))
ggsave(file.path(main_directory, '3D_GAP_vs_NMR_state_2_scatterplot.pdf'), height = 2, width = 2)


# plot intrinsic hydrolysis vs state 2
int_NMR_spread <- data %>% 
  select(-sd) %>% 
  filter(measure %in% c('int', 'NMR')) %>% 
  spread(measure, value) 
int_NMR_spread <- data %>% 
  select(mutant, measure, sd) %>% 
  filter(measure == 'int') %>% 
  inner_join(., int_NMR_spread, by = 'mutant') %>% 
  filter(! is.na(NMR))
intNMRfit <- lm(int ~ NMR, data = int_NMR_spread)
int_NMR_spread %>% 
  ggplot(aes(y = int, x = NMR)) + 
  geom_point(size = 2, alpha = 0.5) +
  geom_errorbar(aes(ymin = int - sd, ymax = int + sd), width = 2, size = 0.1) + 
  geom_text_repel(aes(label = mutant), size = 2, point.padding = 0.2, segment.size = 0.3) +
  ylab(expression(" relative intrinsic GTP hydrolysis rate [s"^-{}^{1}*']')) +
  xlab(expression("\n% "*gamma*" phosphate in "*gamma*" state 2")) +
  theme_classic() +
  geom_abline(intercept = intNMRfit$coefficients[1], slope = intNMRfit$coefficients[2], color = ucsf_colors$pink1, alpha = 0.5) +  
  theme(text = element_text(size = 6),
        axis.title.x = element_text(size = 6), 
        axis.title.y = element_text(size = 6),
        axis.line = element_line(size = 0.1))
ggsave(file.path(extended_directory, 'Ext_Fig6B_intrinsic_vs_NMR_state_2_scatterplot.pdf'), height = 2.5, width = 2.5)


int_GAP_spread <- data %>% 
  select(-sd) %>% 
  filter(measure %in% c('int', 'kcat_Km')) %>% 
  spread(measure, value)
outliers_int_GAP_spread <- int_GAP_spread %>% 
  filter(mutant %in% c('K132H', 'R108I')) %>% 
  mutate('int' = int, 'kcat_Km' = kcat_Km) %>% 
  summarise('ymin' = min(int) - abs(0.1 * min(int)), 'ymax' = max(int) + abs(0.1 * max(int)), 'xmin' = min(kcat_Km) - abs(0.1 * min(kcat_Km)), 'xmax' = max(kcat_Km) + abs(0.1 * max(kcat_Km)))
filtered_int_GAP_spread <- int_GAP_spread %>% 
  filter(! mutant %in% c('K132H'))
intGAPfit <- lm(int ~ kcat_Km, data = filtered_int_GAP_spread)
int_GAP_spread %>% 
  ggplot(aes(y = int, x = kcat_Km)) + 
  geom_point(size = 2, alpha = 0.5) +
  geom_rect(data = outliers_int_GAP_spread, inherit.aes = F, fill = ucsf_colors$gray3, alpha = 0.2,
            aes(xmin = outliers_int_GAP_spread$xmin, xmax = outliers_int_GAP_spread$xmax, ymin = outliers_int_GAP_spread$ymin, ymax = outliers_int_GAP_spread$ymax)) + 
  geom_text_repel(aes(label = mutant), size = 2, point.padding = 0.2, segment.size = 0.3) +
  ylab(expression("intrinsic GTP hydrolysis ln(rate"['(MUTANT)']*"/ rate"['(WILD TYPE)']*"\n")) +
  xlab(expression("GAP mediated GTP hydrolysis relative k"['cat']*"/K"['m']*" ln(k"['cat']*"/K"['m (MUTANT)']*"/ k"['cat']*"/K"['m (WILD TYPE)'])) +
  theme_classic() +
  ggtitle('GAP mediated GTP hydrolysis versus intrinsic GTP hydrolysis') +
  geom_abline(intercept = intGAPfit$coefficients[1], slope = intGAPfit$coefficients[2], color = ucsf_colors$pink1, alpha = 0.3) +  
  theme(text = element_text(size = 6),
        axis.title.x = element_text(size = 6), 
        axis.title.y = element_text(size = 6),
        axis.line = element_line(size = 0.1))
ggsave(file.path(supplemental_directory, 'intrinsic_vs_GAP_scatterplot.pdf'), height = 4, width = 4)

