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
  #mutate('value' = gamma2/gamma1, 'measure' = 'NMR', 'se' = NA) %>%
  #mutate('value' = value/(filter(., variant == 'WT')$value)) %>% 
  mutate('value' = 100*gamma2/(gamma1+gamma2), 'measure' = 'NMR', 'se' = NA) %>%
  select('mutant' = variant, se, value, measure)

# Read in intrinsic and GAP-mediated GTP hydrolysis data
intrinsic_hydrolysis <- read_tsv('Data/RanGAP_assay/intrinsic_hydrolysis.txt') %>% 
  select(mutant, 'int' = mean_rel_rate, 'se' = se_rel_rate) %>% 
  gather("measure", "value", -mutant, -se) %>% 
  unique()

GAP.effic <- read_tsv("Data/RanGAP_assay/GAP_kinetics_MichaelisMenten_parameters.txt") %>% 
  select(mutant, 'kcat_Km' = mean_kcat_Km, se) %>% 
  gather("measure", "value", -mutant, -se) %>% 
  unique()
GAP.kcat <- read_tsv("Data/RanGAP_assay/GAP_kinetics_MichaelisMenten_parameters.txt") %>% 
  select(mutant, 'kcat' = mean_kcat, 'se' = kcat_se) %>% 
  gather("measure", "value", -mutant, -se) %>% 
  unique()
GAP.Km <- read_tsv("Data/RanGAP_assay/GAP_kinetics_MichaelisMenten_parameters.txt") %>% 
  select(mutant, 'Km' = mean_Km, 'se' = Km_se) %>% 
  gather("measure", "value", -mutant, -se) %>% 
  unique()
data <- bind_rows(GAP.effic, intrinsic_hydrolysis, nmr_data, GAP.Km, GAP.kcat)

WT_parameters <- data %>% 
  filter(mutant == "WT") %>% 
  mutate('value'= ifelse(measure == 'NMR', 1, value)) %>% 
  select(-mutant)

rel_data <- data %>% 
  inner_join(., WT_parameters, by = 'measure') %>% 
  select(mutant, measure, 'value' = value.x, 'se' = se.x, 'wt_value' = value.y, 'wt_se' = se.y) %>% 
  mutate('rel_value' = value/wt_value, 'rel_se' = (sqrt( (se/value)^2 + (wt_se/wt_value)^2 ) * value/wt_value)) %>% 
  select(mutant, measure, 'value' = rel_value, 'se' = rel_se)

# plot GAP vs state 2
GAP_NMR_spread <- rel_data %>%
  select(-se) %>%
  filter(measure %in% c('kcat_Km', 'NMR')) %>%
  spread(measure, value)
GAP_NMR_spread <- rel_data %>% 
  select(mutant, measure, se) %>% 
  filter(measure == 'kcat_Km') %>% 
  inner_join(., GAP_NMR_spread, by = 'mutant') %>% 
  filter(! is.na(NMR))
rel_lin_lin_plot <- GAP_NMR_spread %>% 
  ggplot(aes(y = kcat_Km, x = NMR)) + 
  geom_point(size = 2, alpha = 0.5) +
  geom_errorbar(aes(ymin = kcat_Km - se, ymax = kcat_Km + se), width = 2, size = 0.1) + 
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
rel_lin_lin_plot
#ggsave(filename = str_c(main_directory, 'relative_GAP_vs_NMR.pdf'), height = 2.3, width = 2.2)


GAP_NMR_spread <- rel_data %>%
  select(-se) %>%
  filter(measure %in% c('kcat', 'NMR')) %>%
  spread(measure, value)
GAP_NMR_spread <- rel_data %>% 
  select(mutant, measure, se) %>% 
  filter(measure == 'kcat') %>% 
  inner_join(., GAP_NMR_spread, by = 'mutant') %>% 
  filter(! is.na(NMR))
kcat_plot <- GAP_NMR_spread %>% 
  ggplot(aes(y = kcat, x = NMR)) + 
  geom_point(size = 2, alpha = 0.5) +
  geom_errorbar(aes(ymin = kcat - se, ymax = kcat + se), width = 2, size = 0.1) + 
  geom_text_repel(aes(label = mutant), size = 2, point.padding = 0.1, segment.size = 0.3) +
  ylab(expression("(WT/MUT) - GAP mediated GTP hydrolysis relative k"['cat'])) +
  xlab(expression("\n% "*gamma*" phosphate in "*gamma*" state 2")) +
  theme_classic() +
  ggtitle('GAP versus linear gamma 2') +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  theme(text = element_text(size = 6),
        axis.title = element_text(size = 6), 
        axis.text = element_text(size = 6),
        axis.line = element_line(size = 0.1), 
        axis.ticks = element_line(size = 0.1))
kcat_plot

GAP_NMR_spread <- rel_data %>%
  select(-se) %>%
  filter(measure %in% c('Km', 'NMR')) %>%
  spread(measure, value)
GAP_NMR_spread <- rel_data %>% 
  select(mutant, measure, se) %>% 
  filter(measure == 'Km') %>% 
  inner_join(., GAP_NMR_spread, by = 'mutant') %>% 
  filter(! is.na(NMR))
Km_plot <- GAP_NMR_spread %>% 
  ggplot(aes(y = Km, x = NMR)) + 
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = Km - se, ymax = Km + se), width = 2, size = 0.1) + 
  #geom_text_repel(aes(label = mutant), size = 2, point.padding = 0.1, segment.size = 0.3) +
  theme_classic() +
  #ylim(c(0, 6)) +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  theme(text = element_text(size = 6),
        axis.title = element_text(size = 6), 
        axis.text = element_text(size = 6),
        axis.line = element_line(size = 0.1), 
        axis.ticks = element_line(size = 0.1))
Km_plot

linear_fit <- GAP_NMR_spread %>% 
  filter(! mutant %in% c('K132H', 'R78K', 'D79S')) %>% 
  lm(data = ., Km~NMR)
Km_plot_zoom <- GAP_NMR_spread %>% 
  filter(! mutant %in% c('K132H', 'R78K', 'D79S')) %>% 
  ggplot(aes(y = Km, x = NMR)) + 
  geom_point(size = 2, alpha = 0.8) +
  geom_errorbar(aes(ymin = Km - se, ymax = Km + se), width = 2, size = 0.1) + 
  geom_text_repel(aes(label = mutant), size = 2, point.padding = 0.1, segment.size = 0.3) +
  ylab(expression("(WT/MUT) - GAP mediated GTP hydrolysis relative kM")) +
  xlab(expression("\n% "*gamma*" phosphate in "*gamma*" state 2")) +
  theme_classic() +
  geom_abline(intercept = linear_fit$coefficients[1], slope = linear_fit$coefficients[2], color = ucsf_colors$pink1, alpha = 0.5) +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  theme(text = element_text(size = 6),
        axis.title = element_text(size = 6), 
        axis.text = element_text(size = 6),
        axis.line = element_line(size = 0.1), 
        axis.ticks = element_line(size = 0.1))
Km_plot_zoom
plot_grid(Km_plot_zoom) + draw_plot(Km_plot, x = 0.4, y = 0.5, width = 0.6, height = 0.52)
ggsave(filename = str_c(main_directory, "/Km.pdf"), height = 2, width = 2.2)



# plot intrinsic hydrolysis vs state 2
int_NMR_spread <- rel_data %>% 
  select(-se) %>% 
  filter(measure %in% c('int', 'NMR')) %>% 
  spread(measure, value) 
int_NMR_spread <- rel_data %>% 
  select(mutant, measure, se) %>% 
  filter(measure == 'int') %>% 
  inner_join(., int_NMR_spread, by = 'mutant') %>% 
  filter(! is.na(NMR))
intNMRfit <- lm(int ~ NMR, data = int_NMR_spread)
int_NMR_spread %>% 
  ggplot(aes(y = int, x = NMR)) + 
  geom_point(size = 2, alpha = 0.5) +
  geom_errorbar(aes(ymin = int - se, ymax = int + se), width = 2, size = 0.1) + 
  geom_text_repel(aes(label = mutant), size = 2, point.padding = 0.2, segment.size = 0.3) +
  ylab(expression(" relative intrinsic GTP hydrolysis rate [s"^-{}^{1}*']')) +
  xlab(expression("\n% "*gamma*" phosphate in "*gamma*" state 2")) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 100, by = 20)) +
  #geom_abline(intercept = intNMRfit$coefficients[1], slope = intNMRfit$coefficients[2], color = ucsf_colors$pink1, alpha = 0.5) +  
  theme(text = element_text(size = 6),
        axis.title.x = element_text(size = 6), 
        axis.title.y = element_text(size = 6),
        axis.line = element_line(size = 0.1),
        axis.ticks = element_line(size = 0.1))
ggsave(file.path(extended_directory, 'Ext_Fig7B_intrinsic_vs_NMR_state_2_scatterplot.pdf'), height = 2.5, width = 2.5)


int_GAP_spread <- rel_data %>% 
  select(-se) %>% 
  filter(measure %in% c('int', 'kcat_Km')) %>% 
  spread(measure, value)
int_GAP_spread %>% 
  ggplot(aes(y = int, x = kcat_Km)) + 
  geom_point(size = 2, alpha = 0.5) +
  geom_text_repel(aes(label = mutant), size = 2, point.padding = 0.2, segment.size = 0.3) +
  ylab(expression("intrinsic GTP hydrolysis ln(rate"['(MUTANT)']*"/ rate"['(WILD TYPE)']*"\n")) +
  xlab(expression("GAP mediated GTP hydrolysis relative k"['cat']*"/K"['m']*" ln(k"['cat']*"/K"['m (MUTANT)']*"/ k"['cat']*"/K"['m (WILD TYPE)'])) +
  theme_classic() +
  ggtitle('GAP mediated GTP hydrolysis versus intrinsic GTP hydrolysis') +
  theme(text = element_text(size = 6),
        axis.title.x = element_text(size = 6), 
        axis.title.y = element_text(size = 6),
        axis.line = element_line(size = 0.1))
ggsave(file.path(supplemental_directory, 'intrinsic_vs_GAP_scatterplot.pdf'), height = 4, width = 4)

