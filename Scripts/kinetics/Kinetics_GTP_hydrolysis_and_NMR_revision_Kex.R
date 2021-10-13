library(tidyverse)
library(ggrepel)
library(cowplot)
source('ucsf_colors.R')
main_directory <- 'Revisions/Main Figures/Figure3/'
extended_directory <- 'Revisions/Extended_Figures/'
supplemental_directory <- 'Revisions/Supplementary_Files/Supplementary_File_1'

rel_data <- read_tsv('Data/kinetics_data_relative_to_WT.txt')

# plot ln(MUT_GAP/WT_GAP) vs ln(MUT_Kex/WT_Kex) where Kex is gamma2/gamma1 (measure called gamma_ratio)
GAP_temp <- rel_data %>% 
  filter(measure == 'GAP_kcat_Km') %>% 
  select(mutant, "GAP" = ln_rel_to_WT, 'GAP_se' = ln_se)
temp <- rel_data %>% 
  filter(measure == 'gamma_ratio') %>% 
  select(mutant, 'gamma_ratio' = ln_rel_to_WT) %>% 
  inner_join(., GAP_temp, by = 'mutant')
outliers <- temp %>%
  filter(mutant %in% c('R78K', 'D79S', 'K132H')) %>%
  summarise('ymin' = min(GAP) - abs(0.1 * min(GAP)),
            'ymax' = max(GAP) + abs(0.1 * max(GAP)),
            'xmin' = min(gamma_ratio) - abs(0.1 * min(gamma_ratio)),
            'xmax' = max(gamma_ratio) + abs(0.1 * max(gamma_ratio)))
filtered_temp <- temp %>%
  filter(! mutant %in% c('R78K', 'D79S', 'K132H'))
GAPNMRfit <- lm(GAP ~ gamma_ratio, data = filtered_temp)
ln_rel_GAP_gamma2_plot <- temp %>% 
  ggplot(aes(y = GAP, x = gamma_ratio)) + 
  geom_point(size = 2, alpha = 0.5) +
  geom_errorbar(aes(ymin = GAP - GAP_se, ymax = GAP + GAP_se), width = 0.1, size = 0.1) + 
  geom_text_repel(aes(label = mutant), size = 2, point.padding = 0.1, segment.size = 0.3) +
  ylab(expression("(WT/MUT) - GAP mediated GTP hydrolysis relative k"['cat']*"/K"['m'])) +
  #xlab(expression("\n% "*gamma*" phosphate in "*gamma*" state 2")) +
  xlab(expression("\nln(KexMUT/KexWT) "*gamma*" phosphate in "*gamma*" state 2/"*gamma*"state 1")) +
  #scale_x_continuous(breaks = c(-5, -4, -3, -2, -1, 0, 1, 2)) +
  geom_rect(data = outliers, inherit.aes = F, fill = ucsf_colors$gray3, alpha = 0.2,
           aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) + 
  geom_abline(intercept = GAPNMRfit$coefficients[1], slope = GAPNMRfit$coefficients[2], color = ucsf_colors$pink1, alpha = 0.5) +  
  theme_classic() +
  theme(text = element_text(size = 6),
        axis.title = element_text(size = 6), 
        axis.text = element_text(size = 6),
        axis.line = element_line(size = 0.1), 
        axis.ticks = element_line(size = 0.1))
ln_rel_GAP_gamma2_plot

ggsave(filename = str_c(main_directory, 'Fig3d_GAP_kcatKm_vs_ln_NMR_gamma_Kex.pdf'), height = 2.05, width = 2.22)

### Print for Figure source file
temp %>% write_tsv('Per_Figure_source_files/Fig3D.tsv')

### make a plotting function that does just kcat or just Km
plot_MM_parameter_gamma <- function(measure, input_data) {
  plot<- input_data %>% 
  ggplot(aes(y = value, x = gamma_ratio)) + 
  geom_point(size = 2, alpha = 0.5) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 0.1, size = 0.1) + 
  geom_text_repel(aes(label = mutant), size = 2, point.padding = 0.1, segment.size = 0.3) +
  ylab(str_c("ln(", measure, ")")) +
  xlab(expression("\nln(KexMUT/KexWT) "*gamma*" phosphate in "*gamma*" state 2/"*gamma*"state 1")) +
  theme_classic() +
  theme(text = element_text(size = 6),
        axis.title = element_text(size = 6), 
        axis.text = element_text(size = 6),
        axis.line = element_line(size = 0.1), 
        axis.ticks = element_line(size = 0.1))
  return(plot)
}
### plot kcat
kcat_temp <- rel_data %>% 
  filter(measure == 'GAP_kcat') %>% 
  select(mutant, "value" = ln_rel_to_WT, se)
temp <- rel_data %>% 
  filter(measure == 'gamma_ratio') %>% 
  select(mutant, 'gamma_ratio' = ln_rel_to_WT) %>% 
  inner_join(., kcat_temp, by = 'mutant')

plot_MM_parameter_gamma(measure = 'kcatMUT/kcatMUT', input_data = temp)
#ggsave(filename = str_c(extended_directory, 'GAP_kcat_vs_NMR.pdf'), height = 2.3, width = 3)

### plot Km
Km_temp <- rel_data %>% 
  filter(measure == 'GAP_Km') %>% 
  select(mutant, "value" = ln_rel_to_WT, se)
temp <- rel_data %>% 
  filter(measure == 'gamma_ratio') %>% 
  select(mutant, 'gamma_ratio' = ln_rel_to_WT) %>% 
  inner_join(., Km_temp, by = 'mutant')
Km_plot <- plot_MM_parameter_gamma(measure = 'KmMUT/KmWT', input_data = temp)
Km_plot
#ggsave(filename = str_c(extended_directory, 'GAP_Km_vs_NMR.pdf'), height = 2.3, width = 3)


linear_fit <- temp %>%
  filter(! mutant %in% c('K132H', 'R78K', 'D79S')) %>%
  lm(data = ., value~gamma_ratio)
temp_zoom <- temp %>% 
  filter(! mutant %in% c('K132H', 'R78K', 'D79S')) 
Km_plot_zoom <- plot_MM_parameter_gamma(measure = 'KmMUT/KmWT', input_data = temp_zoom)
Km_plot_zoom <- Km_plot_zoom + geom_abline(slope = linear_fit$coefficients[2], intercept = linear_fit$coefficients[1], color = ucsf_colors$pink1, alpha = 0.5)
plot_grid(Km_plot_zoom) + draw_plot(Km_plot, x = 0.5, y = 0.55, width = 0.5, height = 0.5)
#ggsave(filename = str_c(extended_directory, 'GAP_Km_vs_NMR_with_zoom.pdf'), height = 3.5, width = 3.5)



# plot intrinsic hydrolysis vs state 2
int_temp <- rel_data %>% 
  filter(measure == 'int') %>% 
  select(mutant, "int" = ln_rel_to_WT, 'ln_se' = ln_se)
temp <- rel_data %>% 
  filter(measure == 'gamma_ratio') %>% 
  select(mutant, 'gamma_ratio' = ln_rel_to_WT) %>% 
  inner_join(., int_temp, by = 'mutant')
intNMRfit <- lm(int ~ gamma_ratio, data = temp)
temp %>% 
  ggplot(aes(y = int, x = gamma_ratio)) + 
  geom_point(size = 2, alpha = 0.5) +
  geom_errorbar(aes(ymin = int - ln_se, ymax = int + ln_se), width = 0.1, size = 0.1) + 
  geom_text_repel(aes(label = mutant), size = 2, point.padding = 0.2, segment.size = 0.3) +
  ylab(expression(" relative intrinsic GTP hydrolysis rate [s"^-{}^{1}*']')) +
  xlab(expression("\nln(KexMUT/KexWT) "*gamma*" phosphate in "*gamma*" state 2/"*gamma*"state 1")) +
  theme_classic() +
  scale_y_continuous(breaks = c(-3, -2.6, -2.2, -1.8,-1.4, -1, -0.6, -0.2, 0.2, 0.6)) +
  geom_abline(intercept = intNMRfit$coefficients[1], slope = intNMRfit$coefficients[2], color = ucsf_colors$pink1, alpha = 0.5) +  
  theme(text = element_text(size = 6),
        axis.title.x = element_text(size = 6), 
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6),
        axis.line = element_line(size = 0.1),
        axis.ticks = element_line(size = 0.1))
ggsave(file.path(extended_directory, 'EDF_7/EDFig7_intrinsic_vs_NMR_state_2_scatterplot.pdf'), height = 2.5, width = 2.5)

### Print for Figure source file
temp %>% write_tsv('Per_Figure_source_files/EDF8c.tsv')


# plot intrinsic hydrolysis vs GAP-mediated
int_temp <- rel_data %>% 
  filter(measure == 'int') %>% 
  select(mutant, "int" = ln_rel_to_WT, 'ln_se' = ln_se)
temp <- rel_data %>% 
  filter(measure == 'GAP_kcat_Km') %>% 
  select(mutant, 'GAP_kcat_Km' = ln_rel_to_WT) %>% 
  inner_join(., int_temp, by = 'mutant')
intGAPfit <- lm(GAP_kcat_Km ~ int, data = temp)
temp %>% 
  ggplot(aes(x = int, y = GAP_kcat_Km)) + 
  geom_point(size = 2, alpha = 0.5) +
  geom_errorbar(aes(ymin = GAP_kcat_Km - ln_se, ymax = GAP_kcat_Km + ln_se), width = 0.1, size = 0.1) + 
  geom_text_repel(aes(label = mutant), size = 2, point.padding = 0.2, segment.size = 0.3) +
  xlab(expression(" relative intrinsic GTP hydrolysis rate [s"^-{}^{1}*']')) +
  ylab(expression("ln(WT/MUT) - GAP mediated GTP hydrolysis relative k"['cat']*"/K"['m'])) +
  theme_classic() +
  geom_abline(intercept = 0, slope = 1, color = ucsf_colors$pink1, alpha = 0.5) +  
  theme(text = element_text(size = 6),
        axis.title.x = element_text(size = 6), 
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 6), 
        axis.text.y = element_text(size = 6),
        axis.line = element_line(size = 0.1),
        axis.ticks = element_line(size = 0.1))
ggsave(file.path(extended_directory, 'EDF_7/intrinsic_vs_GAP_scatterplot.pdf'), height = 2.5, width = 2.5)
