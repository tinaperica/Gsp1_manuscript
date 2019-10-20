library(tidyverse)
library(RColorBrewer)
source('ucsf_colors.R')
spectral_colors <- brewer.pal(11, name = 'Spectral')
load('Data/spitzemapko_correlations_and_bonferroni_fdr_all.RData')
load('Data/filtered_v5_correlations.RData')
gene_set_data <- read_tsv('Data/gene_sets.txt')
gene_sets <- gene_set_data %>% pull(gene_set) %>% unique()
gene_sets_to_keep <- gene_sets[c(1, 2, 5, 6, 7, 16, 8, 9, 17, 21)]
gene_sets_annotated <- tibble('gene_set' = gene_sets_to_keep,
                              'color' = c(rep('brown', 2), rep(ucsf_colors$orange1, 4), rep(ucsf_colors$cyan1, 4)),
                              'gene_set_color' = c(rep(spectral_colors[10],2),  
                                                   spectral_colors[1:4], spectral_colors[6:9])
                              )
mutants_ordered <- read_tsv('Data/order_of_mutants.txt', col_names = F)
mutants_annotated <- tibble('mutant' = mutants_ordered$X1,
                            'mutant_color' = c(rep(ucsf_colors$orange1, 8), rep(ucsf_colors$cyan1, 13), 'brown'))
write_tsv(mutants_annotated, 'Figure5_Cell_Network/mutant_nodes.txt')
gene_sets_final <- gene_set_data %>% 
  select(query, gene_name, gene_set) %>% 
  inner_join(., gene_sets_annotated) %>% 
  select('allele' = query, gene_set, 'gene_color' = color, gene_set_color)
alleles <- gene_sets_final %>% pull(allele) %>% unique()
gene_sets_final %>% select(allele, gene_color) %>% write_tsv('Figure5_Cell_Network/allele_nodes.txt')
mut_gene_network <- filtered_correlations %>% 
  filter(grepl('GSP1', query_uniq1) & query_uniq2 %in% alleles) %>% 
  separate(query_uniq1, into = c('temp', 'mutant'), sep = ' - ') %>% 
  select(mutant, 'allele' = query_uniq2, pearson, greater_fdr) %>% 
  inner_join(., gene_sets_final, by = 'allele') %>% 
  inner_join(., mutants_annotated, by = 'mutant') %>% 
  mutate('sig' = 1-greater_fdr) %>% 
  select(mutant, allele, pearson, greater_fdr)
write_tsv(mut_gene_network, 'Figure5_Cell_Network/mut_gene_network.txt')


annot <- mutants_annotated %>% 
  select('node' = mutant, 'color' = mutant_color) %>% 
  mutate('gene_set_color' = spectral_colors[11]) %>% 
  bind_rows(., select(gene_sets_final, 'node' = allele, 'color' = gene_color, gene_set_color)) %>% 
  mutate('node_type' = ifelse(node %in% mutants_ordered$X1, 'mutant', 'allele')) %>% 
  mutate('gene_set_color' = ifelse(node %in% mutants_ordered$X1, color, gene_set_color))
all_network <- correlations %>% 
  mutate('query_uniq1' = ifelse(grepl('GSP1', query_uniq1), substring(query_uniq1, first = 8), query_uniq1)) %>% 
  filter((query_uniq1 %in% mutants_ordered$X1 | query_uniq1 %in% alleles) & query_uniq2 %in% alleles) %>% 
  select('node1' = query_uniq1, 'node2' = query_uniq2, pearson, greater_fdr) %>% 
  select(node1, node2, pearson, greater_fdr)
write_tsv(all_network, 'Figure5_Cell_Network/whole_network.txt')
filtered_whole_network <- all_network %>% 
  filter(greater_fdr < 0.05)
write_tsv(filtered_whole_network, 'Figure5_Cell_Network/filtered_whole_network.txt')

write_tsv(annot, 'Figure5_Cell_Network/node_att.txt')

