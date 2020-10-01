library(tidyverse)
load('Data/spitzemapko_correlations_and_bonferroni_fdr_all.RData')

correlations <- correlations %>% 
  filter(grepl('GSP1', query_uniq1) & ! grepl('GSP1', query_uniq2)) %>% 
  filter(! query_uniq1 == 'GSP1 - T34N') %>% 
  select('mutant' = query_uniq1, 'CellMAP_allele' = query_uniq2, 'yeast_gene' = gene_name2, 'Pearson correlation' = pearson, greater_p_value, 'greater_FDR' = greater_fdr, greater_bonferroni)

write_tsv(correlations, 'Revisions/Supplementary_Files/Supplementary_File_3_correlations.txt')

emap_gene_names <- read_tsv('Data/E-MAP/gsp1_pEMAP_avg_merged_gene_names.txt') %>% 
  gather('array allele', 'E-MAP S-score', -Gene) %>% 
  select('query allele name' = Gene, everything())
emap_orf_names <- read_tsv('Data/E-MAP/gsp1_pEMAP_avg_merged.txt') %>% 
  gather('array allele ORF', 'E-MAP S-score', -Gene) %>% 
  mutate('Gene' = 'YLR293C') %>% 
  mutate('array allele ORF' = ifelse(`array allele ORF` == 'WT HIS3D::KAN - WT HIS3D::KAN', 'WT HIS3D::KAN', `array allele ORF`)) %>% 
  mutate('array allele ORF' = ifelse(`array allele ORF` == 'WT LEU2D::KAN-GFP - WT LEU2D::KAN-GFP', 'WT LEU2D::KAN-GFP', `array allele ORF`)) %>% 
  select('query allele ORF' = Gene, 'array allele ORF')
emap <- bind_cols(emap_gene_names, emap_orf_names) %>% 
  select('query allele name (Gsp1 mutant)' = `query allele name`, `query allele ORF`, `array allele`, `array allele ORF`, `E-MAP S-score`) %>% 
  filter(! `query allele name (Gsp1 mutant)` %in% c('GSP1 - T34N', 'GSP1 - NTER3XFLAG WT', 'GSP1 - CTER3XFLAG WT'))
write_tsv(emap, 'Revisions/Supplementary_Files/Supplementary_File_2_Gsp1_E-MAP.txt')
