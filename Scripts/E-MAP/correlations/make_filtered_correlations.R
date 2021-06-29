#### this code makes the filtered_correlations
### it looks for array genes that have significant genetic interactions with at least one of the Gsp1 mutants
### and then keeps only those spitzemapko queries that have significant interactions with at least 10 of these queries
### the point is to only consider pearson correlations that are based on some signal
## it removes all the weak queries and Gsp1 mutants, because they give artificially high correlations
### This filtered correlations output is the basis for all the analysis in the clustering_by_corr directory

### TLD;DR 
### This code outputs two variations of filtered correlations:
### both keep only the queries that have at least 10 significant interactions (outside of the -3/+3 limits) with array genes that have significant interactions with Gsp1 mutants
### 1) corr_of_corr_with_SGA/filtered_v5_correlations.RData keeps only the queries that have at least 1 significant correlation according to Bonfferoni (n = 476)
### 2) corr_of_corr_with_SGA/filtered_v6_correlations.RData keeps only the queries that have at least 2 significant correlation according to Bonfferoni (n = 278)


library(tidyverse)
library(dendextend)
library(factoextra)
source('ucsf_colors.R')
## spitzemapko
load('Data/spitzemapko_for_corr.rda')
orf_to_gene_index <- read_tsv('Data/orf_gene_GO_sgd_annotation.txt', col_names = F) %>% 
  select('ORF' = X1, 'gene_name' = X2)
sga_index <- read_tsv('Data/spitzemapko_query_orf_index.txt') %>% 
  select('query' = query_allele_name, 'ORF' = query_ORF) %>% 
  inner_join(., orf_to_gene_index, by = 'ORF') %>% 
  select(query, gene_name)
#### remove the flag tagged wt strains
spitzemapko_for_corr <- spitzemapko_for_corr %>% 
  filter(! query_allele_name %in% c('GSP1 - NTER3XFLAG WT', 'GSP1 - CTER3XFLAG WT'))


### how many query genes in spitzemapko (including Gsp1 mutants) 3665
spitzemapko_for_corr %>%
  pull(query_allele_name) %>% unique() %>% length()
### chromatin library array genes - 1127 (used for correlation calculations)
chromatin_library_emap_array_genes <- spitzemapko_for_corr %>%
  filter(grepl('GSP1', query_allele_name)) %>%
  pull(array_ORF) %>% unique() 
chromatin_library_emap_array_genes %>% length()

### out of those 1127, 679 have a significant GI with at least one of the Gsp1 mutants (outside of -3, +3)
chromatin_library_emap_array_genes_with_Gsp1_sig_GI <- spitzemapko_for_corr %>%
  filter(grepl('GSP1', query_allele_name)) %>%
  filter((score < -3 | score > 3)) %>%
  pull(array_ORF) %>% unique()
chromatin_library_emap_array_genes_with_Gsp1_sig_GI %>% length()

#### This is now the main step of filtering the query genes - focus on queries that have at least some overlap with our mutants
### for spitzemapko, count for each query how many significant interactions (-3, 3) it has with array genes
### in chromatin_library_emap_array_genes_with_Gsp1_sig_GI
## filter out all that have fewer than 10
filtered_spitzemapko <- spitzemapko_for_corr %>%
  filter((score < -3 | score > 3) &
           array_ORF %in% chromatin_library_emap_array_genes_with_Gsp1_sig_GI ) %>%
  group_by(query_allele_name) %>%
  mutate('count' = n()) %>%
  ungroup() %>%
  select(query_allele_name, array_ORF, score, count) %>%
  arrange(desc(count)) %>%
  filter(count > 10)
top_queries <- filtered_spitzemapko %>% 
  filter(grepl('GSP1', query_allele_name)) %>% 
  select(query_allele_name, count) %>% unique()
#### the query with the highest number of significant interactions is D79A - it has (-3 < S-score > 3) with 369 array genes
### What are the top non-Gsp1 mutants by significant interactions?
top_queries <- filtered_spitzemapko %>% 
  filter(! grepl('GSP1', query_allele_name)) %>% 
  select(query_allele_name, count) %>% unique()
# > top_queries
# A tibble: 3,386 x 2
# query_allele_name count
# <chr>             <int>
# 1 tim54-5002          317
# 2 anp1                295
# 3 ynl181w-5006        273
# 4 mms21-1             270
# 5 ynl181w-5001        269
# 6 sec26-f856aw860a    245
# 7 cdc50               245
# 8 arp3-31             243
# 9 cog3-2              242
# 10 ric1                241
# # … with 3,376 more rows

spitzemapko_for_corr %>% 
  filter(query_allele_name == 'tim54-5002' | query_allele_name == 'anp1') %>% 
  select(-weight) %>% 
  spread(query_allele_name, score) %>% 
  ggplot(aes(x = `tim54-5002`, y = `anp1`)) + geom_point()
### these top_queries seem like deletions that perturb everything and don't seem to be related with our ran mutants
### a lot of significant interactions but no good correlations

#### filtered_spitzemapko has 22 Gsp1 mutants (after removing FLAG WT N and C)
### and in total (including mutants) 3406 queries
filtered_spitzemapko %>% filter(grepl('GSP1', query_allele_name)) %>%
  pull(query_allele_name) %>% unique()
filtered_spitzemapko %>%
  pull(query_allele_name) %>% unique() %>% length()

### since this removed a lot of Gsp1 mutants
### count again significant array genes with leftover mutants: 676
#### (compared to 676 before - only 3 dropped out)
chromatin_library_emap_array_genes_with_Gsp1_sig_GI <- filtered_spitzemapko %>%
  filter(grepl('GSP1', query_allele_name)) %>%
  filter((score < -3 | score > 3)) %>%
  pull(array_ORF) %>% unique()
chromatin_library_emap_array_genes_with_Gsp1_sig_GI %>% length()
### do the same filtering once again
filtered_spitzemapko <- filtered_spitzemapko %>%
  filter((score < -3 | score > 3) &
           array_ORF %in% chromatin_library_emap_array_genes_with_Gsp1_sig_GI ) %>%
  group_by(query_allele_name) %>%
  mutate('count' = n()) %>%
  ungroup() %>%
  select(query_allele_name, array_ORF, score, count) %>%
  arrange(desc(count)) %>%
  filter(count > 10)

####### status check
filtered_spitzemapko$score %>% length()
filtered_spitzemapko %>% pull(query_allele_name) %>% unique() %>% length()
filtered_spitzemapko %>% filter(grepl('GSP1', query_allele_name)) %>%
  pull(query_allele_name) %>% unique()
### filtered_spitzemapko has 181,898 genetic interaction measurements with scores outside the (-3,3) range 
### for 3403 queries 
####  (which includes 22 Gsp1 point mutants)

### based on this decide which queries to keep
### but return all the genetic interactions for those queries (not just the strong ones outside the (-3,3) range
queries_to_keep <- filtered_spitzemapko %>%
  pull(query_allele_name) %>% unique()
filtered_spitzemapko <- spitzemapko_for_corr %>%
  filter(query_allele_name %in% queries_to_keep)
filtered_spitzemapko$score %>% length()
filtered_spitzemapko %>% pull(query_allele_name) %>% unique() %>% length()
##### filtered_spitzemapko has 3,529,000 genetic interactions for 3403 queries (including 22 Gsp1 point mutants)



correlations %>% 
  ggplot(aes(greater_p_value)) + 
  geom_histogram() +
  facet_wrap(~query_uniq1)

#### do the bonferroni and FDR correction for the correlation p-values
correlations <- read_tsv("Data/gsp1_gene_pair_correlations.txt")
correlations <- correlations %>%
  filter(genes %in% queries_to_keep) %>% 
  group_by(mutants) %>% 
  mutate('fdr' = p.adjust(two_sided_p_value, method = 'fdr'),
         'bonferroni' = p.adjust(two_sided_p_value, method = 'bonferroni'),
         'greater_fdr' = p.adjust(greater_p_value, method = 'fdr'),
         'greater_bonferroni' = p.adjust(greater_p_value, method = 'bonferroni')) %>%
  ungroup() %>% 
  rename('query_uniq1' = mutants, 'query_uniq2' = genes)
correlations <- correlations %>%
  inner_join(., sga_index, by = c('query_uniq1' = 'query')) %>%
  select(everything(), 'gene_name1' = gene_name)
correlations <- correlations %>%
  inner_join(., sga_index, by = c('query_uniq2' = 'query')) %>%
  select(everything(), 'gene_name2' = gene_name)
save(correlations, file = 'Data/spitzemapko_correlations_and_bonferroni_fdr_all.RData')
correlations %>% 
  select('mutant' = query_uniq1, 'CellMAP_allele' = query_uniq2, 'yeast_gene' = gene_name2,	
         'Pearson correlation' = 	pearson, 'greater p-value' = greater_p_value, 'greater FDR' = greater_fdr,
         'greater Bonferroni' = greater_bonferroni) %>% 
  arrange(mutant) %>% 
  write_tsv(., 'Data/Supplementary_File_3.txt')
### this data is provided in Supplementary File 3
#### it contains correlations with all queries from filtered_correlations but for ALL mutants, not just the 22





#### filtered correlations is filtered for alleles and Gsp1 mutants
filtered_correlations <- correlations %>%
  filter(query_uniq1 %in% queries_to_keep & query_uniq2 %in% queries_to_keep)
save(filtered_correlations, file = 'Data/filtered_correlations.RData')

### in the next step of filtering add filtering by positive corr FDR/Bonferroni
queries_with_sig_corr <- filtered_correlations %>%   ## 1179
  filter(grepl('GSP1', query_uniq1) & greater_fdr < 0.05) %>% 
  pull(query_uniq2) %>% unique()
queries_with_sig_corr %>% length()
queries_with_two_sig_corr <- filtered_correlations %>% ## 718
  filter(grepl('GSP1', query_uniq1) & greater_fdr < 0.05) %>% 
  group_by(query_uniq2) %>% 
  summarize('sig_count' = n()) %>% 
  filter(sig_count > 1) %>% 
  pull(query_uniq2) %>% unique()
queries_with_two_sig_corr %>% length()
queries_with_bonf_sig_corr <- filtered_correlations %>% ### 380
  filter(grepl('GSP1', query_uniq1) & greater_bonferroni < 0.05) %>% 
  pull(query_uniq2) %>% unique()
queries_with_bonf_sig_corr %>% length()
queries_with_two_sig_bonf <- filtered_correlations %>% ## 224
  filter(grepl('GSP1', query_uniq1) & greater_bonferroni < 0.05) %>% 
  group_by(query_uniq2) %>% 
  summarize('sig_count' = n()) %>% 
  filter(sig_count > 1) %>% 
  pull(query_uniq2) %>% unique()

filtered_correlations_BNF_1 <- filtered_correlations %>% 
  filter(query_uniq2 %in% queries_with_bonf_sig_corr) %>% 
  arrange(query_uniq2)
filtered_correlations_BNF_1 %>% pull(query_uniq2) %>% unique()
#### final filtered_correlations version BNF_1 has 380 queries
save(filtered_correlations_BNF_1, file = 'Data/filtered_correlations_BNF_1.RData')




filtered_correlations_BNF_2 <- filtered_correlations %>% 
  filter(query_uniq2 %in% queries_with_two_sig_bonf)
filtered_correlations_BNF_2 %>% pull(query_uniq2) %>% unique()
#### final filtered_correlations version BNF_2 has 224 queries
save(filtered_correlations_BNF_2, file = 'Data/filtered_correlations_BNF_2.RData')

### filtering still keeps it to 22 Gsp1 mutants
filtered_correlations_BNF_2 %>% filter(grepl('GSP1', query_uniq1) & (! grepl('GSP1', query_uniq2))) %>% pull(query_uniq1) %>% unique()

### for all the downstream analysis use the version BNF_1 or version BNF_2 - 
#### it's a combination of filtered_correlations (only queries with at least 10 sig interactions with arrays that have sig interactions with at least one Ran mutant)
#### and then additional filtering, that the query needs to have at least one significant correlation with any of the mutants (significance based on Boferroni corrected p-value)
