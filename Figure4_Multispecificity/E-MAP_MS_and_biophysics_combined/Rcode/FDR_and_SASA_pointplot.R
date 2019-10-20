library(tidyverse)
source('ucsf_colors.R')
interfaces <- read_tsv('Data/SASA_interfaces.txt')
core_partners <- interfaces %>% pull(partner) %>% unique()
load('Data/spitzemapko_correlations_and_bonferroni_fdr_all.RData')
load('Data/filtered_v4_correlations.RData')
strong_mutants <- filtered_correlations %>% 
  filter(grepl('GSP1', query_uniq1)) %>% 
  pull(query_uniq1) %>% unique()
correlations <- correlations %>% 
  filter(query_uniq1 %in% strong_mutants) %>% 
  select('mutant' = query_uniq1, 'query' = query_uniq2, 'gene_name' = gene_name2, 'FDR' = greater_fdr, pearson) %>% 
  filter(grepl('GSP1', mutant) & gene_name %in% core_partners) %>% 
  mutate('mutant' = substring(mutant, first = 8)) %>% 
  mutate('residue' = as.numeric(substring(mutant, 2, (nchar(mutant)-1))))
queries <- correlations$gene_name %>% unique()
ordered_mutants <- correlations %>% 
  arrange(residue) %>% pull(mutant) %>% unique()
data <- correlations %>% 
  inner_join(., interfaces, by = c('residue' = 'yeastresnum')) %>% 
  select(-interface) %>% 
  complete(partner, nesting(mutant, query, gene_name, FDR, pearson, residue), fill = list(deltarASA = 0))



partners_data <- data %>% 
  filter(gene_name %in% queries & gene_name == partner) %>% 
  select(mutant, gene_name, pearson, partner, deltarASA) %>% 
  mutate('pearson' = ifelse(pearson < -0.3, -0.3, pearson)) %>% 
  mutate('pearson' = ifelse(pearson > 0.3, 0.3, pearson)) %>% 
  mutate('mutant' = factor(mutant, ordered_mutants)) %>% 
  arrange(mutant)

drASA_plot <- partners_data %>% 
  ggplot(aes(x = gene_name, y = mutant, fill = pearson, size = deltarASA)) +
  geom_point(shape = 21, stroke = 0.1) +
  scale_fill_gradientn(colors = c(ucsf_colors$pink1, 'white', ucsf_colors$green1), 
                       limits = c(-0.3, 0.3), breaks=c(-0.3, 0, 0.3), na.value = ucsf_colors$gray3) +
  scale_size(range = c(0, 4.5), breaks = c(0, 0.1, 0.25, 0.5, 1)) +
  labs(size = expression(Delta*'rASA')) +
  theme_classic() +
  ylab('Point mutation in Ran') +
  xlab('\nRan partner') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = 0.05),
        axis.ticks.length = unit(0.05, 'cm'),
        legend.position = 'bottom', 
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(.15, "cm"),
        axis.line = element_line(size = 0.1),
        plot.margin = unit(c(0, 0, 0, 0), 'cm'))

partners_data <- data %>% 
  filter(gene_name %in% core_partners) %>% 
  select(mutant, gene_name, pearson, FDR) %>% 
  unique() %>% 
  mutate('pearson' = ifelse(pearson < -0.3, -0.3, pearson)) %>% 
  mutate('pearson' = ifelse(pearson > 0.3, 0.3, pearson)) %>% 
  mutate('mutant' = factor(mutant, levels = ordered_mutants)) %>% 
  arrange(mutant)


pvalue_plot <- partners_data %>% 
  ggplot(aes(x = gene_name, y = mutant, fill = pearson, size = FDR)) +
  geom_point(shape = 21, stroke = 0.1) +
  scale_fill_gradientn(colors = c(ucsf_colors$pink1, 'white', ucsf_colors$green1), 
                       limits = c(-0.3, 0.3), breaks=c(-0.3, 0, 0.3), na.value = ucsf_colors$gray3) +
  scale_size("FDR", range = c(4.5, 0.1), breaks = c(0, 0.001, 0.05, 0.1, 0.5)) +
  theme_classic() +
  ylab(element_blank()) +
  xlab('Ran partner') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = 0.05),
        axis.ticks.length = unit(0.05, 'cm'),
        legend.position = 'bottom', 
        legend.direction = 'horizontal',
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"),
        axis.line = element_line(size = 0.1),
        plot.margin = unit(c(0, 0, 0, 0), 'cm')) +
  guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5),
         size = guide_legend(title.position = "top", title.hjust = 0.5))


combined_plots <- plot_grid(drASA_plot, pvalue_plot, align = 'h', rel_widths = c(1.1, 1))
pdf(file = 'Figure4/E-MAP_MS_and_biophysics_combined/FDR_and_deltaSASA_plot_partners.pdf', width = 3.5, height = 3.5)
combined_plots
dev.off()
