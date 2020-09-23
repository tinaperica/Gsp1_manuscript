library(tidyverse)
source('ucsf_colors.R')
data <- read_tsv('Data/APMS_data.txt') %>% 
  filter(norm == 'eqM') %>% 
  filter(! Prey_gene_name == 'GSP1') %>% 
  mutate('Prey_gene_name' = str_c(substr(Prey_gene_name, 1, 1), tolower(substr(Prey_gene_name, 2, nchar(Prey_gene_name))))) %>% 
  mutate('sample' = ifelse(tag == 'N', str_c('N-3xFL-', mutant), str_c(mutant, '-C-3xFL'))) %>% 
  mutate('Prey_gene_name' = case_when(is.na(Prey_gene_name) ~ PreyORF,
                                      !is.na(Prey_gene_name) ~ Prey_gene_name)) %>% 
  select(sample, tag, mutant, residue, Prey_gene_name, log2FC, adj.pvalue) %>% 
  unique()
max_value <- data %>% filter(! abs(log2FC) == Inf) %>% pull(log2FC) %>% max(na.rm = T)
min_value <- data %>% filter(! abs(log2FC) == Inf) %>% pull(log2FC) %>% min(na.rm = T)
### set Inf and -Inf values to max and min values for calculating IQR
data <- data %>% 
  mutate('log2FC' = ifelse(log2FC == Inf, max_value, log2FC)) %>% 
  mutate('log2FC' = ifelse(log2FC == -Inf, min_value, log2FC))
interquartile_range <- data %>% 
  group_by(Prey_gene_name) %>% 
  summarize('IQR' = IQR(log2FC, na.rm = T)) %>% 
  arrange(desc(IQR))
interquartile_range_partners <- interquartile_range %>% 
  filter(Prey_gene_name %in% c('Rna1', 'Srm1', 'Yrb1', 'Kap95', 'Srp1', 'Pse1', 'Spa2', 'Pup2')) %>% 
  bind_rows(., tibble('Prey_gene_name' = c('mean', '+1 s.d.'), 
                      'IQR' = c(mean(interquartile_range$IQR), mean(interquartile_range$IQR) + sd(interquartile_range$IQR))))
interquartile_range %>% ggplot(aes(x = IQR)) + 
  geom_histogram(color = ucsf_colors$gray3, fill = ucsf_colors$gray3, alpha = 0.2) +
  theme_classic() +
  labs(color = 'Gsp1 partner as AP-MS prey protein') +
  xlab(expression('Interquartile range of log'[2]*'(abundance'[MUT]*'/abundance'[WT]*') for each prey protein')) +
  geom_segment(data = interquartile_range_partners, 
               aes(x = IQR, xend = IQR, y = 0, yend = 20, color = Prey_gene_name),
               arrow = arrow(length = unit(0.2,"cm"), ends = 'first'),
               size = 1.1, alpha = 0.8) +
  scale_color_manual(limits = c("Rna1", 'Srp1', "Kap95", "Srm1", 'Yrb1', 'Pse1', 'mean', '+1 s.d.'),
                     values=c(ucsf_colors$orange1, ucsf_colors$pink2, ucsf_colors$green1, ucsf_colors$cyan1, ucsf_colors$yellow1, ucsf_colors$green3, ucsf_colors$gray1, ucsf_colors$gray2)) +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.ticks = element_line(size = 0.05),
        axis.ticks.length = unit(0.05, 'cm'),
        legend.position = 'bottom',
        legend.direction = 'horizontal', 
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        #legend.key.size = unit(0.3, "cm"),
        axis.line = element_line(size = 0.1))
ggsave('Revisions/Extended_Figures/Ext_Fig5_IQR_of_log2FC.pdf', width = 4, height = 3) 
