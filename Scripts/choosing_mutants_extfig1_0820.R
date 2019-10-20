library(tidyverse)
source('ucsf_colors.R')
residue_data <- read_tsv('Data/interface_residues_and_mutations.txt')

# prepare partner colors to match extended data figure 1A
partner_colors <- c(ucsf_colors$cyan1, ucsf_colors$orange1, ucsf_colors$purple1,
                    rep(ucsf_colors$blue1, 2), rep(ucsf_colors$yellow1, 2),
                    ucsf_colors$pink1, rep(ucsf_colors$green1, 8))
names(partner_colors) <- c('Srm1','Rna1','Ntf2','Nup1','Nup60','Yrb1','Yrb2','Srp1',
                           'Crm1','Los1','Msn5','Pse1','Kap95','Cse1','Mtr10','Kap104')

# CRM1 has residue 113 listed as core AND rim?
# residue_data %>% 
#   group_by(partner, yeastresnum) %>% 
#   summarise(n = n()) %>% 
#   filter(n > 1)

residue_data %>% 
  filter(nucleotide != 'other') %>% 
  pull(yeastresnum) %>% 
  unique()

residue_data %>% 
  mutate('partner' = str_c(substr(partner, 1, 1), tolower(substr(partner, 2, nchar(partner))))) %>% 
  mutate('yeastresnum' = factor(yeastresnum, unique(residue_data$yeastresnum))) %>% 
  mutate('interface' = factor(interface, c('core', 'rim', 'support', 'none'))) %>% 
  ggplot(aes(x = yeastresnum, y = partner)) + 
  geom_point(aes(size = interface, shape = mutated, color = partner)) +
  scale_colour_manual(values = partner_colors, guide=FALSE) +
  scale_shape_manual(values = c('mutated' = 16, 'other' = 1),
                     labels = c('mutated', 'not mutated'),) +
  scale_size_manual(values = c('core' = 2.5, 'rim' = 1.5, 'support' = 1.5, 'none' = 0.1)) +
  xlab("Gsp1 interface residues") + ylab("Gsp1 interaction partner") +
  theme_bw() +
  theme(
    text = element_text(family = "Helvetica", size = 6),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
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
    labs(shape = '', color = '')

ggsave('Extended_Figures/Ext_Fig1B_Overview_of_residue_selection_for_Gsp1_PM_0820.pdf', width = 7.2, height = 3)
dev.off()

