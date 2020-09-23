library(tidyverse)
library(cowplot)
source('ucsf_colors.R')
data <- read_tsv('Data/RanGAP_assay/simulated_data_error_estimate/simulated_fits_errors.txt')

### kcat plot
ordered_samples_by_kcat <- data %>% arrange(exp_kcat) %>% pull(Data)
kcat_plot <- data %>% 
  mutate('Data' = factor(Data, ordered_samples_by_kcat)) %>% 
  arrange(Data) %>% 
  ggplot(aes(x = Data, y = exp_kcat)) +
  ylab('kcat') +
  geom_point(aes(x = Data, y = fit_kcat), color = ucsf_colors$navy1, alpha = 0.75) +
  geom_point(color = ucsf_colors$orange1, alpha = 0.75) +
  geom_errorbar(aes(ymin = fit_kcat - RMSD_kcat, ymax = fit_kcat + RMSD_kcat), 
                width = 0.6, size = 0.3, alpha = 0.9, color = ucsf_colors$navy2) +
  theme_classic() +
  theme(text = element_text(size = 8, family='Helvetica'), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 8), 
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.key.size =  unit(0.5, 'cm'),
        legend.position = 'rigth',
        axis.line = element_line(size = 0.1),
        axis.ticks = element_line(size = 0.1)
        )
ggsave('~/Desktop/kcat_exp_vs_fit.pdf', width = 7, height = 3)
### Km
ordered_samples_by_Km <- data %>% arrange(exp_Km) %>% pull(Data)
Km_plot <- data %>% 
  mutate('Data' = factor(Data, ordered_samples_by_Km)) %>% 
  arrange(Data) %>% 
  ggplot(aes(x = Data, y = exp_Km)) +
  ylab('Km') +
  geom_point(aes(x = Data, y = fit_Km), color = ucsf_colors$navy1, alpha = 0.75, size = 0.5) +
  geom_point(color = ucsf_colors$orange1, alpha = 0.75, size = 0.5) +
  geom_errorbar(aes(ymin = fit_Km - RMSD_Km, ymax = fit_Km + RMSD_Km), 
                width = 0.6, size = 0.3, alpha = 0.9, color = ucsf_colors$navy2) +
  theme_classic() +
  theme(text = element_text(size = 6, family='Helvetica'), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 6), 
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size =  unit(0.22, 'cm'),
        legend.position = 'none',
        axis.line = element_line(size = 0.1),
        axis.ticks = element_line(size = 0.1)
  )

# kcat/Km
ordered_samples_by_kcat_Km <- data %>% arrange(exp_kcat_Km) %>% pull(Data)
kcat_Km_plot <- data %>% 
  mutate('Data' = factor(Data, ordered_samples_by_kcat_Km)) %>% 
  arrange(Data) %>% 
  ggplot(aes(x = Data, y = exp_kcat_Km)) +
  ylab('kcat_Km') +
  geom_point(aes(x = Data, y = fit_kcat_Km), color = ucsf_colors$navy1, alpha = 0.75, size = 0.5) +
  geom_point(color = ucsf_colors$orange1, alpha = 0.75, size = 0.5) +
  geom_errorbar(aes(ymin = fit_kcat_Km - RMSD_kcat_Km, ymax = fit_kcat_Km + RMSD_kcat_Km), 
                width = 0.6, size = 0.3, alpha = 0.9, color = ucsf_colors$navy2) +
  theme_classic() +
  theme(text = element_text(size = 6, family='Helvetica'), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 6), 
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size =  unit(0.22, 'cm'),
        legend.position = 'none',
        axis.line = element_line(size = 0.1),
        axis.ticks = element_line(size = 0.1)
  )

## error (%RMSD) versus C (C = [S]/Km)
percent_RMSD_plot <- data %>% 
  select(`C = [S]/Km`, 'kcat' = percent_deviation_kcat, 'Km' = percent_deviation_Km, 'kcat/Km' = percent_deviation_kcat_Km) %>% 
  pivot_longer(., col = c(kcat, Km, `kcat/Km`),
               names_to = '% RMSD', values_to = 'value') %>%
  ggplot(aes(x = `C = [S]/Km`, y = value, color = `% RMSD`)) +
  geom_point(alpha = .75) +
  scale_color_manual(values = c(ucsf_colors$green1, ucsf_colors$yellow1, ucsf_colors$blue1)) +
  ylab('% RMSD') +
  xlab('[S]/Km') +
  theme_classic() +
  theme(text = element_text(size = 6, family = 'Helvetica'), 
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6), 
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size =  unit(0.22, 'cm'),
        legend.position = 'right',
        axis.line = element_line(size = 0.1),
        axis.ticks = element_line(size = 0.1)
  )
  
#### percent standard deviation
percent_sd_plot <- data %>% 
  select(`C = [S]/Km`, 'kcat' = fit_kcat_sd_percent, 'Km' = fit_Km_sd_percent, 'kcat/Km' = fit_kcat_Km_sd_percent) %>% 
  pivot_longer(., col = c(kcat, Km, `kcat/Km`),
               names_to = '% Std.Dev', values_to = 'value') %>%
  ggplot(aes(x = `C = [S]/Km`, y = value, color = `% Std.Dev`)) +
  geom_point(alpha = 0.75) +
  scale_color_manual(values = c(ucsf_colors$green1, ucsf_colors$yellow1, ucsf_colors$blue1)) +
  ylab('% Std.Dev') +
  xlab('[S]/Km') +
  theme_classic() +
  theme(text = element_text(size = 6, family = 'Helvetica'), 
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6), 
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size =  unit(0.22, 'cm'),
        legend.position = 'none',
        axis.line = element_line(size = 0.1),
        axis.ticks = element_line(size = 0.1)
  )
plot_grid(percent_RMSD_plot, percent_sd_plot, align = 'v', nrow = 1)






rmsd <- data %>% 
  select(`C = [S]/Km`, 'kcat' = percent_deviation_kcat, 'Km' = percent_deviation_Km, 'kcat/Km' = percent_deviation_kcat_Km,) %>% 
  pivot_longer(., col = c(kcat, Km, `kcat/Km`),
               names_to = 'deviation from mean', values_to = 'value') %>% 
  mutate('measure' = 'RMSD')
sd <- data %>% 
  select(`C = [S]/Km`, 'kcat' = fit_kcat_sd_percent, 'Km' = fit_Km_sd_percent, 'kcat/Km' = fit_kcat_Km_sd_percent) %>% 
  pivot_longer(., col = c(kcat, Km, `kcat/Km`),
               names_to = 'deviation from mean', values_to = 'value') %>% 
  mutate('measure' = 'Std.Dev')
deviation_data <- bind_rows(rmsd, sd)

deviation_data %>%
  ggplot(aes(x = `C = [S]/Km`, y = value, color = `deviation from mean`)) +
  geom_point(alpha = .75) +
  facet_wrap(~measure) +
  scale_color_manual(values = c(ucsf_colors$green1, ucsf_colors$yellow1, ucsf_colors$blue1)) +
  ylab('% deviation from the mean (RMSD or standard deviaton)') +
  xlab('[S]/Km') +
  theme_classic() +
  theme(text = element_text(size = 7, family = 'Helvetica'), 
        axis.text.x = element_text(size = 7), 
        axis.text.y = element_text(size = 7), 
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size =  unit(0.3, 'cm'),
        legend.position = 'right',
        axis.line = element_line(size = 0.1),
        axis.ticks = element_line(size = 0.1)
  )
ggsave('Revisions/Supplementary_Files/Supplementary_File_1/Supplementary_File_1_Figures/Sup_Fig_14_GAP_IMM_error_estimation.pdf', width = 7.5, height = 3.5)

deviation_data %>% group_by(`deviation from mean`, measure) %>% 
  summarise('max' = max(value))
