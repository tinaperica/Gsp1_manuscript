#### plot point plot with fold change
library(tidyverse)
library(cowplot)
library(extrafont)
loadfonts()
source('ucsf_colors.R')

data <- read_tsv('Data/APMS_data.txt') %>% 
  filter(norm == 'eqM') %>% 
  select(tag, mutant, residue, Prey_gene_name, log2FC, adj.pvalue, deltarASA, interface_partner) %>% 
  unique() %>% 
  complete(interface_partner, nesting(tag, mutant, residue, Prey_gene_name, log2FC, adj.pvalue), fill = list(deltarASA = 0))
# core_partners
core_partners <- data %>% pull(interface_partner) %>% unique()
ordered_mutants <- data %>% select(residue, mutant) %>% unique() %>% arrange(residue) %>% pull(mutant)
### preys
preys <- data %>% pull(Prey_gene_name) %>% unique()


partners_data <- data %>% 
  filter(Prey_gene_name %in% core_partners & interface_partner %in% preys & Prey_gene_name == interface_partner) %>% 
  select(mutant, tag, Prey_gene_name, log2FC, interface_partner, deltarASA) %>% 
  mutate('sample' = ifelse(tag == 'N', str_c('N-3xFL-', mutant), str_c(mutant, '-C-3xFL'))) %>% 
  mutate('log2FC' = ifelse(log2FC < -4, -4, log2FC)) %>% 
  mutate('log2FC' = ifelse(log2FC > 4, 4, log2FC)) %>% 
  mutate('mutant' = factor(mutant, levels = ordered_mutants)) %>% 
  arrange(mutant)
partners_N_data <- partners_data %>% 
  filter(tag == 'N') %>% 
  mutate('mutant' = factor(mutant, levels = ordered_mutants)) %>% 
  arrange(mutant)
partners_C_data <- partners_data %>% 
  filter(tag == 'C') %>% 
  mutate('mutant' = factor(mutant, levels = ordered_mutants)) %>% 
  arrange(mutant)

drASA_plot <- partners_data %>% 
  ggplot(aes(x = Prey_gene_name, y = mutant, color = log2FC, size = deltarASA, shape = tag)) +
  geom_point() +
  scale_color_gradientn(limits = c(-4, 4), breaks = c(-4, 0, 4), 
                        colors = c(ucsf_colors$pink1, 'white', ucsf_colors$blue1),
                        guide = 'none') +
  scale_shape_manual(values = c("\u25D7", "\u25D6"), guide = 'none') +
  scale_size(range = c(0, 6), breaks = c(0, 0.25, 0.5, 0.75)) +
  labs(size = expression(Delta*'rASA')) +
  theme_classic() +
  ylab('Gsp1 point mutant') +
  xlab('\nGsp1 partner by AP-MS') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = 0.05),
        axis.ticks.length = unit(0.05, 'cm'),
        legend.position = 'bottom', 
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        #legend.key.size = unit(.15, "cm"),
        axis.line = element_line(size = 0.1),
        plot.margin = unit(c(0, 0, 0, 0), 'cm'))

partners_data <- data %>% 
  filter(Prey_gene_name %in% core_partners & interface_partner %in% preys) %>% 
  select(mutant, tag, Prey_gene_name, log2FC, adj.pvalue) %>% 
  unique() %>% 
  mutate('sample' = ifelse(tag == 'N', str_c('N-3xFL-', mutant), str_c(mutant, '-C-3xFL'))) %>% 
  mutate('log2FC' = ifelse(log2FC < -4, -4, log2FC)) %>% 
  mutate('log2FC' = ifelse(log2FC > 4, 4, log2FC)) %>% 
  mutate('mutant' = factor(mutant, levels = ordered_mutants)) %>% 
  arrange(mutant)
partners_N_data <- partners_data %>% 
  filter(tag == 'N') %>% 
  mutate('mutant' = factor(mutant, levels = ordered_mutants)) %>% 
  arrange(mutant)
partners_C_data <- partners_data %>% 
  filter(tag == 'C') %>% 
  mutate('mutant' = factor(mutant, levels = ordered_mutants)) %>% 
  arrange(mutant)

pvalue_plot <- partners_data %>% 
  ggplot(aes(x = Prey_gene_name, y = mutant, color = log2FC, size = adj.pvalue, shape = tag)) +
  geom_point() +
  scale_color_gradientn(limits = c(-4, 4), breaks = c(-4, 0, 4), 
                        colors = c(ucsf_colors$pink1, 'white', ucsf_colors$blue1)) +
  scale_shape_manual(values = c("\u25D7", "\u25D6"), guide = 'none') +
  labs(color = expression('log'[2]*'(abundance'['MUT']*'/abundance'['WT']*')')) +
  scale_size("FDR", range = c(5, 0), breaks = c(0, 0.001, 0.05, 0.1)) +
  theme_classic() +
  ylab(element_blank()) +
  xlab('\nGsp1 partner by AP-MS') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = 0.05),
        axis.ticks.length = unit(0.05, 'cm'),
        legend.position = 'bottom', 
        legend.direction = 'horizontal',
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        #legend.key.size = unit(0.3, "cm"),
        axis.line = element_line(size = 0.1),
        plot.margin = unit(c(0, 0, 0, 0), 'cm')) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
          size = guide_legend(title.position="top", title.hjust = 0.5))


combined_plots <- plot_grid(drASA_plot, pvalue_plot, align = 'h', rel_widths = c(1.38, 1))
cairo_pdf(file = 'Figure2_APMS/AP-MS_semicircle_plot_partners.pdf', family = "Arial Unicode MS", width = 2.2, height = 3.8)
combined_plots
dev.off()

cairo_pdf(file = 'Figure2_APMS/AP-MS_semicircle_plot_partners_for_legend.pdf', family = "Arial Unicode MS", width = 6, height = 4)
combined_plots
dev.off()


gene_sets <- read_tsv('Data/gene_sets.txt') %>% 
  select(gene_name, gene_set, Description) %>% unique()
data_gene_set <- data %>% 
  select(-interface_partner, -deltarASA) %>%
  unique() %>% 
  inner_join(., gene_sets, by = c('Prey_gene_name' = 'gene_name')) %>% 
  arrange(gene_set)
data_gene_set %>% pull(Prey_gene_name) %>% unique()
