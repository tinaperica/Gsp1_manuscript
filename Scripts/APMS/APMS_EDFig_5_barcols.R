library(tidyverse)
source('ucsf_colors.R')

interfaces <- read_tsv('Data/Gsp1_interfaces_SASA_and_conservation.txt') %>% 
  filter(protein == 'GSP1') %>% 
  mutate('interface_partner' = str_c(substr(partner, 1, 1), tolower(substr(partner, 2, nchar(partner))))) %>% 
  select(interface_partner, interface, deltarASA, deltaASA, rASAc, 'residue' = yeast_num)

apms_data <- read_tsv('Data/APMS_data.txt') %>% 
  filter(! Prey_gene_name == 'GSP1') %>% 
  filter(norm == 'eqM') %>% 
  mutate('Prey_gene_name' = str_c(substr(Prey_gene_name, 1, 1), tolower(substr(Prey_gene_name, 2, nchar(Prey_gene_name))))) %>% 
  mutate('interface_partner' = str_c(substr(interface_partner, 1, 1), tolower(substr(interface_partner, 2, nchar(interface_partner))))) %>% 
  mutate('Prey_gene_name' = case_when(is.na(Prey_gene_name) ~ PreyORF,
                                      !is.na(Prey_gene_name) ~ Prey_gene_name)) %>% 
  mutate('sample' = ifelse(tag == 'N', str_c('N-3xFL-', mutant), str_c(mutant, '-C-3xFL'))) %>% 
  select(sample, tag, mutant, residue, Prey_gene_name, log2FC, adj.pvalue, interface_partner) %>% 
  unique() %>% 
  complete(interface_partner, nesting(sample, tag, mutant, residue, Prey_gene_name, log2FC, adj.pvalue), fill = list(deltarASA = 0))

data <- apms_data %>% 
  inner_join(interfaces, by = c('interface_partner', 'residue')) %>% 
  arrange(sample, interface_partner, residue, Prey_gene_name) %>% 
  complete(interface_partner, nesting(sample, tag, mutant, residue, Prey_gene_name, log2FC, adj.pvalue), fill = list(deltarASA = 0)) %>% 
  mutate('interface' = ifelse((interface == 'rim' | interface == 'support'), 'rim/support', interface)) %>% 
  mutate('interface' = ifelse(is.na(interface), 'not_interface', interface)) 


order_of_samples <- data %>% arrange(residue) %>% pull(sample) %>% unique()
order_of_preys <- c('Rna1', 'Srm1', 'Yrb1', 'Srp1', 'Kap95', 'Pse1')
data %>% 
  mutate('sample' = factor(sample, order_of_samples)) %>% 
  mutate('Prey_gene_name' = factor(Prey_gene_name, order_of_preys)) %>% 
  arrange(sample, Prey_gene_name) %>% 
  filter(interface_partner == Prey_gene_name & interface == 'core') %>% 
  unique() %>% 
  ggplot(aes(x = sample, y = log2FC, fill = Prey_gene_name)) + 
  geom_col(color = 'black', size = 0.25) +
  scale_fill_manual(name = 'AP-MS Prey proteins when mutation is\nin the core with the prey interface', values = c(ucsf_colors$orange1, ucsf_colors$cyan1, ucsf_colors$yellow1, ucsf_colors$pink1, ucsf_colors$green1, ucsf_colors$green3)) +
  theme_classic() +
  xlab('\namino or carboxy terminally FLAG-tagged point mutants of Gsp1 (AP-MS Bait)') +
  theme(
    text = element_text(family = "Helvetica", size = 6),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 90,
                               # color =ifelse(grepl('N', order_of_samples), "red", "blue") # for coloring x-axis labels
                               ),
    axis.ticks = element_line(size = 0.05),
    axis.ticks.length = unit(0.05, 'cm'),
    legend.position = 'bottom',
    legend.direction = 'horizontal',
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.key.size = unit(0.35, "cm"),
    axis.line = element_line(size = 0.1)
  )
ggsave('Final_Formatting/temp/EDF_APMS_barcols.pdf', height = 2.9, width = 5.5)
### print for source file
data %>% write_tsv('Per_Figure_source_files/EDF6_AB.txt')
#####

data %>% 
  mutate('sample' = factor(sample, order_of_samples)) %>% 
  mutate('Prey_gene_name' = factor(Prey_gene_name, order_of_preys)) %>% 
  arrange(sample, Prey_gene_name) %>% 
  filter(interface_partner == Prey_gene_name) %>% 
  unique() %>% 
  ggplot(aes(x = sample, y = log2FC, fill = Prey_gene_name)) + 
  geom_col(color = 'black', size = 0.25) +
  scale_fill_manual(name = 'AP-MS Prey proteins when mutation is\nin the core with the prey interface', values = c(ucsf_colors$orange1, ucsf_colors$cyan1, ucsf_colors$yellow1, ucsf_colors$pink1, ucsf_colors$green1, ucsf_colors$green3)) +
  theme_classic() +
  xlab('\namino or carboxy terminally FLAG-tagged point mutants of Gsp1 (AP-MS Bait)') +
  theme(
    text = element_text(family = "Helvetica", size = 6),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    axis.text.x = element_text(angle = 90),
    axis.ticks = element_line(size = 0.05),
    axis.ticks.length = unit(0.05, 'cm'),
    legend.position = 'bottom',
    legend.direction = 'horizontal',
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.key.size = unit(0.35, "cm"),
    axis.line = element_line(size = 0.1)
  )
ggsave('Final_Formatting/temp/EDF_APMS_barcols_all.pdf', height = 2.9, width = 5.5)

data %>% 
  mutate('sample' = factor(sample, order_of_samples)) %>% 
  mutate('Prey_gene_name' = factor(Prey_gene_name, order_of_preys)) %>% 
  arrange(sample, Prey_gene_name) %>% 
  filter(Prey_gene_name == 'Srm1' & sample == 'N-3xFL-T34A')

