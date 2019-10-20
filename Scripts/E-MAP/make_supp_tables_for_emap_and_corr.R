load('Data/spitzemapko_correlations_and_bonferroni_fdr_all.RData')

correlations <- correlations %>% 
  filter(grepl('GSP1', query_uniq1) & ! grepl('GSP1', query_uniq2)) %>% 
  select('mutant' = query_uniq1, 'CellMAP_allele' = query_uniq2, 'yeast_gene' = gene_name2, 'Pearson correlation' = pearson, greater_p_value, 'greater_FDR' = greater_fdr, greater_bonferroni)

write_tsv(correlations, 'Supplementary_Data_Tables/Supp_Table4_correlations.txt')

emap <- read_tsv('Data/E-MAP/gsp1_pEMAP_avg_merged_gene_names.txt') %>% 
  gather('library mutant', 'E-MAP S-score', -Gene) %>% 
  select('Gsp1 point mutant' = Gene, everything()) %>% 
  arrange(`Gsp1 point mutant`, `library mutant`)
write_tsv(emap, 'Supplementary_Data_Tables/Supp_Table3_E-MAP_S-scores.txt')
