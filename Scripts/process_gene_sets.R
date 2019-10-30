### this script combines the gene_sets.txt file and the 7 clusters from Fig. 4b
### to make the supplementary table with gene sets information
library(tidyverse)
gene_sets <- read_tsv('Data/gene_sets.txt')
seven_clusters <- read_csv('Supplementary_Data_Tables/Excel_files/corr_clustering_heatmap_cluster_gene_sets.csv')
final_table <- gene_sets %>% 
  select(-cluster) %>% 
  left_join(., seven_clusters, by = c('query' = 'strain')) %>% 
  select(query, cluster, gene_name, gene_set) %>% 
  mutate('cluster' = ifelse(is.na(cluster), 'expanded dataset', cluster)) %>% 
  select(query, gene_name, gene_set, cluster)

write_tsv(final_table, path = 'Supplementary_Data_Tables/Excel_files/gene_sets_final.txt')
