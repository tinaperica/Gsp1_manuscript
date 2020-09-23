library(tidyverse)
library(ggrepel)
library(cowplot)
source('ucsf_colors.R')
main_directory <- 'Revisions/Main Figures/Figure3/'
extended_directory <- 'Revisions/Extended_Figures/'
supplemental_directory <- 'Revisions/Supplementary_Files/Supplementary_File_1'

# Read in NMR data
nmr_data <- read_csv('Data/31P_NMR_Data.csv', col_types = cols()) %>%
  filter(peak %in% c('gamma1','gamma2')) %>%
  select(variant, peak, integral_rel) %>%
  spread(peak, integral_rel) %>%
  mutate('gamma1' = case_when(is.na(gamma1) ~ 0.0, T ~ gamma1),
         'gamma2' = case_when(is.na(gamma2) ~ 0.0, T ~ gamma2)) %>%
  mutate('value' = gamma2/gamma1, 'measure' = 'NMR', 'sd' = NA) %>%
  mutate('value' = value/(filter(., variant == 'WT')$value)) %>% 
  # mutate('value' = 100*gamma2/(gamma1+gamma2), 'measure' = 'NMR', 'sd' = NA) %>%
  select('mutant' = variant, sd, value, measure)

  
# Read in intrinsic and GAP-mediated GTP hydrolysis data
intrinsic_hydrolysis <- read_tsv('Data/RanGAP_assay/intrinsic_hydrolysis.txt') %>% 
  select(mutant, 'int'= mean_log_rel_rate, 'sd' = sd_log_rel_rate) %>% 
  gather("measure", "value", -mutant, -sd) %>% 
  unique()

MM.data <- read_tsv("Data/RanGAP_assay/GAP_kinetics_MichaelisMenten_parameters.txt") %>% 
  select(mutant, 'kcat_Km' = mean_kcat_Km, sd) %>% 
  gather("measure", "value", -mutant, -sd) %>% 
  unique()
data <- bind_rows(MM.data, intrinsic_hydrolysis, nmr_data)

WT_parameters <- data %>% 
  filter(mutant == "WT") %>% 
  mutate('value'= ifelse(measure == 'NMR', 1, value)) %>% 
  select(-mutant)

rel_data <- data %>% 
  inner_join(., WT_parameters, by = 'measure') %>% 
  select(mutant, measure, 'value' = value.x, 'sd' = sd.x, 'wt_value' = value.y, 'wt_sd' = sd.y) %>% 
  mutate('rel_value' = value/wt_value,
         'rel_sd' = (sqrt( (sd/value)^2 + (wt_sd/wt_value)^2 ) * value/wt_value)) %>% 
  select(mutant, measure, 'value' = rel_value, 'sd' = rel_sd)

# plot GAP vs state 2
GAP_NMR_spread <- rel_data %>%
  select(-sd) %>%
  filter(measure %in% c('kcat_Km', 'NMR')) %>%
  spread(measure, value)
GAP_NMR_spread <- rel_data %>% 
  select(mutant, measure, sd) %>% 
  filter(measure == 'kcat_Km') %>% 
  inner_join(., GAP_NMR_spread, by = 'mutant') %>% 
  filter(! is.na(NMR))

# plot

# log(rel. GAP efficiency) vs. log(rel. K eq)
data <- GAP_NMR_spread %>% 
  filter(!is.infinite(NMR)) %>% 
  filter(NMR != 0)

data %>% 
  ggplot(aes(y = log(kcat_Km), x = log(NMR))) +
  geom_smooth(method='lm', formula= y~x, se=F,
              data = subset(GAP_NMR_spread, ! mutant %in% c('R78K', 'D79S'))) +
  geom_point(size = 2, alpha = 0.5) +
  # geom_errorbar(aes(ymin = kcat_Km - sd, ymax = kcat_Km + sd), width = 2, size = 0.1) + 
  geom_text_repel(aes(label = mutant), size = 2, point.padding = 0.1, segment.size = 0.3) +
  ylab(expression("\n(WT/MUT) - GAP mediated GTP hydrolysis relative k"['cat']*"/K"['m'])) +
  xlab(expression(" log(K_eq MUT) - log(K_eq WT)")) +
  theme_classic() +
  ggtitle('log(GAP) versus log(K)') +
  # scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  theme(text = element_text(size = 6),
        axis.title = element_text(size = 6), 
        axis.text = element_text(size = 6),
        axis.line = element_line(size = 0.1), 
        axis.ticks = element_line(size = 0.1))
ggsave(filename = str_c(main_directory, 'log_relative_GAP_vs_NMR_logK.pdf'), height = 2.3, width = 2)


# rel. GAP efficiency) vs. rel. K eq

data <- filter(GAP_NMR_spread, !is.infinite(NMR))
data %>% 
  ggplot(aes(y = kcat_Km, x = NMR)) + 
  geom_smooth(method='lm', formula= y~x, se=F,
              data = subset(GAP_NMR_spread, ! mutant %in% c('R78K', 'D79S'))) +
  geom_point(size = 2, alpha = 0.5) +
  # geom_errorbar(aes(ymin = kcat_Km - sd, ymax = kcat_Km + sd), width = 2, size = 0.1) + 
  geom_text_repel(aes(label = mutant), size = 2, point.padding = 0.1, segment.size = 0.3) +
  ylab(expression("(WT/MUT) - GAP mediated GTP hydrolysis relative k"['cat']*"/K"['m'])) +
  xlab(expression("K_eq (MUT) / K_eq (WT)")) +
  theme_classic() +
  ggtitle('linear GAP versus linear K') +
  # scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  theme(text = element_text(size = 6),
        axis.title = element_text(size = 6), 
        axis.text = element_text(size = 6),
        axis.line = element_line(size = 0.1), 
        axis.ticks = element_line(size = 0.1))

ggsave(filename = str_c(main_directory, 'lin_relative_GAP_vs_NMR_K.pdf'), height = 2.3, width = 2)

# Km vs. K plot
# Read in intrinsic and GAP-mediated GTP hydrolysis data
intrinsic_hydrolysis <- read_tsv('Data/RanGAP_assay/intrinsic_hydrolysis.txt') %>% 
  select(mutant, 'int'= mean_log_rel_rate, 'sd' = sd_log_rel_rate) %>% 
  gather("measure", "value", -mutant, -sd) %>% 
  unique()
MM.data <- read_tsv("Data/RanGAP_assay/GAP_kinetics_MichaelisMenten_parameters.txt") %>% 
  select(mutant, "Km" = mean_Km, sd) %>% 
  gather("measure", "value", -mutant, -sd) %>% 
  unique()
data <- bind_rows(MM.data, intrinsic_hydrolysis, nmr_data)

WT_parameters <- data %>% 
  filter(mutant == "WT") %>% 
  mutate('value'= ifelse(measure == 'NMR', 1, value)) %>% 
  select(-mutant)

rel_data <- data %>% 
  inner_join(., WT_parameters, by = 'measure') %>% 
  select(mutant, measure, 'value' = value.x, 'sd' = sd.x, 'wt_value' = value.y, 'wt_sd' = sd.y) %>% 
  mutate('rel_value' = value/wt_value,
         'rel_sd' = (sqrt( (sd/value)^2 + (wt_sd/wt_value)^2 ) * value/wt_value)) %>% 
  select(mutant, measure, 'value' = rel_value, 'sd' = rel_sd)

GAP_NMR_spread <- rel_data %>%
  select(-sd) %>%
  filter(measure %in% c('Km', 'NMR')) %>%
  spread(measure, value)
GAP_NMR_spread <- rel_data %>% 
  select(mutant, measure, sd) %>% 
  filter(measure == 'Km') %>% 
  inner_join(., GAP_NMR_spread, by = 'mutant') %>% 
  filter(! is.na(NMR))


# GAP Km vs. log(rel. K eq)
data <- GAP_NMR_spread %>% 
  filter(!is.infinite(NMR)) %>% 
  filter(NMR != 0)

data %>% 
  ggplot(aes(y = Km, x = log(NMR))) + 
  geom_smooth(method='lm', formula= y~x, se=F,
              data = subset(GAP_NMR_spread, ! mutant %in% c('R78K', 'D79S'))) +
  geom_point(size = 2, alpha = 0.5) +
  # geom_errorbar(aes(ymin = kcat_Km - sd, ymax = kcat_Km + sd), width = 2, size = 0.1) + 
  geom_text_repel(aes(label = mutant), size = 2, point.padding = 0.1, segment.size = 0.3) +
  ylab('GAP Km') +
  xlab(expression("log(K_eq (MUT) / K_eq (WT))")) +
  theme_classic() +
  ggtitle('GAP Km versus log(K)') +
  # scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  theme(text = element_text(size = 6),
        axis.title = element_text(size = 6), 
        axis.text = element_text(size = 6),
        axis.line = element_line(size = 0.1), 
        axis.ticks = element_line(size = 0.1))

ggsave(filename = str_c(main_directory, 'GAP_Km_vs_NMR_logK.pdf'), height = 2.3, width = 2)

# log(GAP Km) vs. log(rel. K eq)

data %>% 
  ggplot(aes(y = log(Km), x = log(NMR))) + 
  geom_smooth(method='lm', formula= y~x, se=F,
              data = subset(GAP_NMR_spread, ! mutant %in% c('R78K', 'D79S'))) +
  geom_point(size = 2, alpha = 0.5) +
  # geom_errorbar(aes(ymin = kcat_Km - sd, ymax = kcat_Km + sd), width = 2, size = 0.1) + 
  geom_text_repel(aes(label = mutant), size = 2, point.padding = 0.1, segment.size = 0.3) +
  ylab('GAP Km') +
  xlab(expression("log(K_eq (MUT) / K_eq (WT))")) +
  theme_classic() +
  ggtitle('GAP log(Km) versus log(K)') +
  # scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  theme(text = element_text(size = 6),
        axis.title = element_text(size = 6), 
        axis.text = element_text(size = 6),
        axis.line = element_line(size = 0.1), 
        axis.ticks = element_line(size = 0.1))

ggsave(filename = str_c(main_directory, 'log_GAP_Km_vs_NMR_logK.pdf'), height = 2.3, width = 2)


# # compute correlation coefficient
# # NOTE CANNOT COMPUTE EXACT P-VALUE BECAUSE OF TIES
# # including outliers
# spearman <-
#   cor.test(GAP_NMR_spread$NMR, GAP_NMR_spread$kcat_Km, method='spearman') %>% 
#   `$`(estimate) %>% 
#   round(3)
# 
# # without outliers
# spearman_no_outliers <-
#   cor.test(filter(GAP_NMR_spread, ! mutant %in% c('R78K', 'D79S', 'K132H'))$NMR,
#            filter(GAP_NMR_spread, ! mutant %in% c('R78K', 'D79S', 'K132H'))$kcat_Km,
#            method='spearman') %>% 
#   `$`(estimate) %>% 
#   round(3)
# 
# spearman_text = paste0(
#   'Spearman rank corr: \n    ', spearman,
#   '\nSpearman rank corr w/out outliers: \n    ', spearman_no_outliers
# )
# 
# # annotate the plot with the spearman coefficients
# corr_annotation <- annotate(
#   'text', label = spearman_text, x = 0, y = 2.75, hjust = 0, family='Helvetica', size = 2)
#                             

# rel_lin_lin_plot + corr_annotation

rel_lin_lin_plot
ggsave(filename = str_c(main_directory, 'relative_GAP_vs_NMR_K.pdf'), height = 2.3, width = 2)



GAP_NMR_spread <- data %>%
  select(-sd) %>%
  filter(measure %in% c('kcat_Km', 'NMR')) %>%
  spread(measure, value)
GAP_NMR_spread <- data %>% 
  select(mutant, measure, sd) %>% 
  filter(measure == 'kcat_Km') %>% 
  inner_join(., GAP_NMR_spread, by = 'mutant') %>% 
  filter(! is.na(NMR))
lin_lin_plot <- GAP_NMR_spread %>% 
  ggplot(aes(y = kcat_Km, x = NMR)) + 
  geom_point(size = 2, alpha = 0.5) +
  geom_errorbar(aes(ymin = kcat_Km - sd, ymax = kcat_Km + sd), width = 2, size = 0.1) + 
  geom_text_repel(aes(label = mutant), size = 2, point.padding = 0.1, segment.size = 0.3) +
  ylab(expression("(WT/MUT) - GAP mediated GTP hydrolysis relative k"['cat']*"/K"['m'])) +
  xlab(expression("\n% "*gamma*" phosphate in "*gamma*" state 2")) +
  theme_classic() +
  ggtitle('GAP versus linear gamma 2') +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  theme(text = element_text(size = 6),
        axis.title = element_text(size = 6), 
        axis.text = element_text(size = 6),
        axis.line = element_line(size = 0.1), 
        axis.ticks = element_line(size = 0.1))
lin_lin_plot
ggsave(filename = str_c(main_directory, 'absolute_GAP_vs_NMR.pdf'), height = 2, width = 2)
