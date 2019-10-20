library(tidyverse)
source('ucsf_colors.R')
residue_data <- read_tsv('Data/interface_residues_and_mutations.txt')
residue_data %>% 
  mutate('partner' = str_c(substr(partner, 1, 1), tolower(substr(partner, 2, nchar(partner))))) %>% 
  mutate('yeastresnum' = factor(yeastresnum, unique(residue_data$yeastresnum))) %>% 
  mutate('interface' = factor(interface, c('core', 'rim', 'support', 'none'))) %>% 
  ggplot(aes(x = yeastresnum, y = partner)) + 
  geom_point(aes(size = interface, shape = nucleotide, color = mutated, alpha = switch)) +
  scale_colour_manual(breaks = 'mutated', values = c(ucsf_colors$navy1,  ucsf_colors$gray3)) +
  scale_shape_manual(breaks = 'within 5 A of the nucleotide', values = c(16, 17)) +
  scale_alpha_manual(breaks = "switch I or switch II", values = c(1, 0.4)) +
  scale_size_manual(values = c('core' = 3.5, 'rim' = 2, 'support' = 2, 'none' = 0.1)) +
  xlab("Gsp1 interface residues") + ylab("Gsp1 interaction partner") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = 0.05),
        axis.ticks.length = unit(0.05, 'cm'),
        axis.line = element_line(size = 0.1),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.spacing = unit(0.01, "cm"),
        legend.margin = margin(0, 0.5, 0, 0.5, unit = "cm"),
        legend.position = 'bottom') +
  labs(shape = '', alpha = '', color = '')

ggsave('Extended_Figures/Ext_Fig1B_Overview_of_residue_selection_for_Gsp1_PM.pdf', width = 7.2, height = 3)
dev.off()

