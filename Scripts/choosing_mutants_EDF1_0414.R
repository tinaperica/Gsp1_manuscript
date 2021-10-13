library(tidyverse)
source('ucsf_colors.R')
interfaces <- read_tsv('Data/Gsp1_interfaces_SASA_and_conservation.txt') %>% 
  filter(protein == 'GSP1') %>% 
  select(partner, interface, yeast_num)
nucleotide_interactions <- read_tsv('Data/interface_residues_and_mutations.txt') %>% 
  select(partner, 'yeast_num' = yeastresnum, nucleotide, mutated, switch)

residue_data <- interfaces %>% 
  full_join(., nucleotide_interactions, by = c('partner', 'yeast_num')) %>% 
  mutate('interface' = ifelse(interface == 'not_interface', 'none', interface)) %>% 
  mutate('nucleotide' = ifelse(is.na(nucleotide), 'other', nucleotide)) %>% 
  mutate('mutated' = ifelse(is.na(mutated), 'other', mutated)) %>% 
  mutate('switch' = ifelse(is.na(switch), 'other', switch)) %>% 
  arrange(yeast_num, partner)
residues_to_keep_in_plot <- residue_data %>% 
  filter(interface == 'core' | mutated == 'mutated' | nucleotide != 'other' | switch != 'other') %>% 
  pull(yeast_num) %>% unique()
residue_data <- residue_data %>% 
  filter(yeast_num %in% residues_to_keep_in_plot) %>% 
  complete(., yeast_num, nesting(partner)) %>% 
  arrange(yeast_num, partner) %>% 
  mutate('interface' = ifelse(is.na(interface), 'none', interface)) %>% 
  mutate('nucleotide' = ifelse(is.na(nucleotide), 'other', nucleotide)) %>% 
  mutate('mutated' = ifelse(is.na(mutated), 'other', mutated)) %>% 
  mutate('switch' = ifelse(is.na(switch), 'other', switch))
  

# prepare partner colors to match extended data figure 1A
partner_colors <- c(ucsf_colors$cyan1, ucsf_colors$orange1, ucsf_colors$purple1,
                    rep(ucsf_colors$blue1, 2), rep(ucsf_colors$yellow1, 2),
                    ucsf_colors$pink1, rep(ucsf_colors$green1, 8))
names(partner_colors) <- c('Srm1','Rna1','Ntf2','Nup1','Nup60','Yrb1','Yrb2','Srp1',
                           'Crm1','Los1','Msn5','Pse1','Kap95','Cse1','Mtr10','Kap104')

residue_data %>% 
  filter(nucleotide != 'other') %>% 
  pull(yeast_num) %>% 
  unique()

residue_data %>% 
  mutate('partner' = str_c(substr(partner, 1, 1), tolower(substr(partner, 2, nchar(partner))))) %>% 
  mutate('yeast_num' = factor(yeast_num, unique(residue_data$yeast_num))) %>% 
  mutate('interface' = factor(interface, c('core', 'rim', 'support', 'none'))) %>% 
  arrange(yeast_num) %>% 
  ggplot(aes(x = yeast_num, y = partner)) + 
  geom_point(aes(size = interface, shape = mutated, color = partner)) +
  scale_colour_manual(values = partner_colors, guide = FALSE) +
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

ggsave('Revisions/Extended_Figures/Ext_Fig1B_Overview_of_residue_selection_for_Gsp1_PM_20200414.pdf', width = 7.4, height = 3)


# save for Per Figure source file
residue_data %>% 
  mutate('partner' = str_c(substr(partner, 1, 1), tolower(substr(partner, 2, nchar(partner))))) %>% 
  mutate('yeast_num' = factor(yeast_num, unique(residue_data$yeast_num))) %>% 
  mutate('interface' = factor(interface, c('core', 'rim', 'support', 'none'))) %>% 
  arrange(yeast_num) %>% 
  write_tsv('Per_Figure_source_files/EDF1g.txt')
