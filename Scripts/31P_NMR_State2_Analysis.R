library(tidyverse)
library("ggpubr")
source('ucsf_colors.R')
# set ordering of variants to be plotted
ordering = c('H141R','K132H','Y157A','D79S','WT',
             'G80A','R108L','R78K','T34G','Q147E',
             'T34L','T34A','T34Q','T34E')

data <- read_csv('Data/31P_NMR/31P_NMR_Data.csv', col_types = cols()) %>%
  filter(peak %in% c('gamma1','gamma2')) %>%
  select(variant, peak, integral_rel) %>%
  spread(peak, integral_rel) %>%
  mutate('gamma1' = case_when(is.na(gamma1) ~ 0.0, T ~ gamma1),
         'gamma2' = case_when(is.na(gamma2) ~ 0.0, T ~ gamma2)) %>%
  mutate('percent_state_2' = 100*gamma2/(gamma1+gamma2)) %>%
  select(variant, percent_state_2) %>%
  mutate(variant = factor(variant, levels = ordering))

# to plot the barplot of state 2 values
data %>% 
  ggplot(aes(variant, percent_state_2,
             fill = factor(ifelse(variant=="WT","Wild Type", "Mutant")))) +
  geom_bar(color = "black", stat = "identity") +
  coord_flip(clip = "off") +
  scale_fill_manual(name = "variant", values = c(ucsf_colors$gray2, ucsf_colors$navy2)) +
  geom_text(aes(label = round(percent_state_2, 1), hjust = -0.1), size=7) +
  xlab('') +
  ylab('Percent in State 2') +
  ylim(0,100) +
  theme(text = element_text(size = 28, color='black'),
        axis.text.y = element_text(color = 'black', size = 20),
        axis.text.x = element_text(color = 'black', size = 20),
        axis.title.x = element_text(size = 25),
        legend.position = 'none',
        plot.margin = unit(c(0, 1.5, 0, 0),'cm'),
        panel.background = element_rect(fill = 'transparent'))

ggsave("Biophysics_Figure3/NMR_plots/horizontal_nmr_barplot.pdf", height = 8, width = 4, bg = 'transparent')


data %>% 
  ggplot(aes(variant, percent_state_2,
             fill = factor(ifelse(variant=="WT","Wild Type", "Mutant")))) +
  geom_bar(color = "black", stat = "identity") +
  coord_flip(clip = "off") +
  scale_fill_manual(name = "variant", values = c(ucsf_colors$gray2, ucsf_colors$navy2)) +
  geom_text(aes(label = round(percent_state_2, 1), hjust = -0.1), size=7) +
  xlab('') +
  ylab('Percent in State 2') +
  ylim(0,100) +
  theme(text = element_text(size = 28, color='black'),
        axis.text.y = element_text(color = 'black', size = 20),
        axis.text.x = element_text(color = 'black', size = 20),
        axis.title.x = element_text(size = 25),
        legend.position = 'none',
        #plot.margin = unit(c(0.5,1.5,0.5,0.5),'cm'),
        plot.margin = unit(c(0, 1.5, 0, 0),'cm'),
        panel.background = element_rect(fill = 'transparent'))

ggsave("Biophysics_Figure3/NMR_plots/vertical_nmr_barplot.pdf", height = 8, width = 4, bg = 'transparent')
