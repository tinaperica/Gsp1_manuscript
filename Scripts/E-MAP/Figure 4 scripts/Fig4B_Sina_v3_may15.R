
##### load libraries
library(tidyverse)
library(ggforce)
source('ucsf_colors.R')

##### load datafiles
# load('Data/spitzemapko_correlations_and_bonferroni_fdr_all.RData')
load('Data/filtered_v6_correlations.RData')
name2ORF <- read_delim('Data/spitzemap_name2ORF_index.txt', delim = '\t', col_types = cols())

# set mutant groups based on clustering
mutant_group_order = c('I', 'II', 'III')
mutant_groups <- 
  list('I' = c('D79S', 'T34Q', 'T34E', 'K101R', 'T34G', 'D79A'),
       'II' = c('T34A', 'Q147E', 'R108I', 'R108L', 'G80A', 'Y157A', 'H141E'),
       'III' = c('H141R', 'R108Y', 'R108Q', 'R108G', 'Y148I', 'H141I', 'R112A', 'R112S', 'R78K')) %>% 
  stack() %>% 
  rename('mutant' = 'values', 'mutant_group' = 'ind') %>% 
  mutate(mutant_group = factor(mutant_group, levels = mutant_group_order))

# clean correlations dataset, so each row is a correlation between a mutant and a strain
corr_for_sina <-
  filtered_correlations %>%
  filter(grepl('GSP1', query_uniq1), !grepl('GSP1', query_uniq2)) %>% 
  separate(query_uniq1, sep = ' - ', into = c('GSP1','query_uniq1')) %>%
  filter(query_uniq1 %in% unlist(mutant_groups)) %>% 
  select(query_uniq1, query_uniq2, pearson, greater_fdr) %>%
  rename('mutant' = query_uniq1, 'strain' = query_uniq2) %>% 
  left_join(mutant_groups, by = 'mutant') 


# annotate strains based on gene set  
gene_sets_order <- c('nuclear pore complex', 'spindle assembly checkpoint', 'tRNA modification', 'other')
gene_sets_colors <- c(ucsf_colors$green1, ucsf_colors$pink1, ucsf_colors$blue1, ucsf_colors$gray3)

gene_sets <- 
  'Supplementary_Data_Tables/Excel_files/gene_sets_final.txt' %>% 
  read_delim(delim = '\t', col_types = cols()) %>% 
  select(query, gene_set) %>% 
  rename('strain' = query) %>% 
  filter(gene_set %in% gene_sets_order)

data <-
  corr_for_sina %>% 
  left_join(gene_sets, by = 'strain') %>% 
  filter(greater_fdr < 0.05) %>% 
  mutate(gene_set = ifelse(is.na(gene_set), 'other', gene_set) %>% 
           factor(levels = gene_sets_order))

ggplot(data, aes(x = mutant_group, y = pearson, size = greater_fdr,)) +
  geom_violin(data = filter(data, gene_set == 'other'),
              fill = ucsf_colors$gray3, color=NA,
              width = 0.8, alpha = 0.2) +
  geom_jitter(data = filter(data, gene_set == 'nuclear pore complex'),
              mapping = aes(x = as.numeric(mutant_group) - 0.3, size = greater_fdr, color = gene_set),
              width = 0.1) +
  geom_jitter(data = filter(data, gene_set == 'spindle assembly checkpoint'),
              mapping = aes(size = greater_fdr, color = gene_set),
              width = 0.1) +
  geom_jitter(data = filter(data, gene_set == 'tRNA modification'),
              mapping = aes(x = as.numeric(mutant_group) + 0.3, size = greater_fdr, color = gene_set),
              width = 0.1) +
  
  scale_color_manual(name='Gene set', values=gene_sets_colors) +
  scale_size_continuous(name = 'P-value', range = c(0.5, 0.0),
                        breaks = c(0.001, 0.01, 0.03, 0.05), limits = c(0, 0.06)) +
  ylim(c(0.05, 0.45)) + xlab('Gsp1 mutant group') +
  ylab('Gsp1 mutant - gene\nPearson correlation') +
  guides(size=guide_legend(nrow=4,byrow=TRUE,
                           override.aes = list(fill=NA)),
         color=guide_legend(nrow=3,byrow=TRUE)) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica", size = 6),
    axis.title = element_text(size = 6), axis.text = element_text(size = 6),
    axis.ticks = element_line(size = 0.05), axis.ticks.length = unit(0.05, 'cm'),
    # legend.position = 'top',
    legend.position = 'none',
    legend.spacing.y = unit(0.01, 'cm'),
    legend.box = 'horizontal', legend.box.just = 'left',
    legend.text = element_text(size = 6), legend.title = element_text(size = 6),
    legend.margin = margin(t = 0, unit='cm'),
    plot.margin = margin(t = 0, unit='cm'),
    axis.line = element_line(size = 0.1),
    strip.text.x = element_text(size = 6)
  )

ggsave('Revisions/Main Figures/Figure4/Fig4_Sina.pdf', height = 1.5, width = 2.5)
dev.off()

