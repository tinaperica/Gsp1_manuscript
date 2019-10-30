#### this code makes the filtered_correlations
### it looks for array genes that have significant genetic itneractions with at least one of the Gsp1 mutants
### and then keeps only those spitzemapko queries that have siginificant interactions with at least 10 of these queries
### the point is to only consider pearson correlations that are based on some signal
## it removes all the weak queries and Gsp1 mutants, because they give arificially high correlations
### This filtered correlations output is the basis for all the analysis in the clustering_by_corr directory

### TLD;DR 
### This code outputs two variations of filtered correlations:
### both keep only the queries that have at least 10 siginificant interactions (outside of the -3/+3 limits) with array genes that have significant interactions with Gsp1 mutants
### 1) corr_of_corr_with_SGA/filtered_v5_correlations.RData keeps only the queries that have at least 1 siginificant correlation according to Bonfferoni (n = 476)
### 2) corr_of_corr_with_SGA/filtered_v6_correlations.RData keeps only the queries that have at least 1 siginificant correlation according to Bonfferoni (n = 278)
library(tidyverse)
library(dendextend)
library(factoextra)
source('ucsf_colors.R')
## spitzemapko
load('basic_E-MAP_data/spitzemapko_for_corr.rda')
orf_to_gene_index <- read_tsv('orf_gene_GO_sgd_annotation.txt', col_names = F) %>% 
  select('ORF' = X1, 'gene_name' = X2)
sga_index <- read_tsv('basic_E-MAP_data/spitzemapko_query_orf_index.txt') %>% 
  select('query' = query_allele_name, 'ORF' = query_ORF) %>% 
  inner_join(., orf_to_gene_index, by = 'ORF') %>% 
  select(query, gene_name)
#### remove the flag tagged wt strains
spitzemapko_for_corr <- spitzemapko_for_corr %>% 
  filter(! query_allele_name %in% c('GSP1 - NTER3XFLAG WT', 'GSP1 - CTER3XFLAG WT'))

#### do the bonferroni and FDR correction  for the correlation p-values - RUN THIS ONCE - uncomment this sectio  to run
# load('corr_of_corr_with_SGA/correlations/spitzemapko_correlations_and_pvalues_all.RData')
# correlations <- correlations %>%
#   group_by(query_uniq1) %>%
#   mutate('fdr' = p.adjust(two_sided_p_value, method = 'fdr'),
#          'bonferroni' = p.adjust(two_sided_p_value, method = 'bonferroni'),
#          'greater_fdr' = p.adjust(greater_p_value, method = 'fdr'),
#          'greater_bonferroni' = p.adjust(greater_p_value, method = 'bonferroni')) %>% 
#   ungroup()
# correlations <- correlations %>%
#   inner_join(., sga_index, by = c('query_uniq1' = 'query')) %>%
#   select(everything(), 'gene_name1' = gene_name)
# correlations <- correlations %>% 
#   inner_join(., sga_index, by = c('query_uniq2' = 'query')) %>%
#   select(everything(), 'gene_name2' = gene_name)
# save(correlations, file = 'corr_of_corr_with_SGA/correlations/spitzemapko_correlations_and_bonferroni_fdr_all.RData')
load('corr_of_corr_with_SGA/correlations/spitzemapko_correlations_and_bonferroni_fdr_all.RData')  ### this data is provided in Supplementary Data Table 5
### how many query genes in spitzemapko (including Gsp1 mutants) 3665
spitzemapko_for_corr %>%
  pull(query_allele_name) %>% unique() %>% length()
### chromatin library array genes - 1129 (used for correlation calculations)
chromatin_library_emap_array_genes <- spitzemapko_for_corr %>%
  filter(grepl('GSP1', query_allele_name)) %>%
  pull(array_ORF) %>% unique()

### out of those 1129, 680 have a significant GI with at least one of the Gsp1 mutants (outside of -3, +3)
chromatin_library_emap_array_genes_with_Gsp1_sig_GI <- spitzemapko_for_corr %>%
  filter(grepl('GSP1', query_allele_name)) %>%
  filter((score < -3 | score > 3)) %>%
  pull(array_ORF) %>% unique()

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
#   1 tim54-5002          317
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
### and in total (including mutants) 3408 queries
filtered_spitzemapko %>% filter(grepl('GSP1', query_allele_name)) %>%
  pull(query_allele_name) %>% unique()
filtered_spitzemapko %>%
  pull(query_allele_name) %>% unique() %>% length()

### since this removed a lot of Gsp1 mutants
### count again significant array genes with leftover mutants: 677 
#### (compared to 680 before - only 3 dropped out)
chromatin_library_emap_array_genes_with_Gsp1_sig_GI <- filtered_spitzemapko %>%
  filter(grepl('GSP1', query_allele_name)) %>%
  filter((score < -3 | score > 3)) %>%
  pull(array_ORF) %>% unique()

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

#### double filtered_spitzemapko now has 22 Gsp1 mutant
### and in total (including mutants) 3405 queries (16 queries dropped out)
filtered_spitzemapko %>% filter(grepl('GSP1', query_allele_name)) %>%
  pull(query_allele_name) %>% unique()
queries_to_keep <- filtered_spitzemapko %>%
  pull(query_allele_name) %>% unique()

filtered_spitzemapko <- spitzemapko_for_corr %>%
  filter(query_allele_name %in% queries_to_keep)

#### filtered correlations are all queries that have at least 10 array gene overlaps
filtered_correlations <- correlations %>%
  filter(query_uniq1 %in% queries_to_keep & query_uniq2 %in% queries_to_keep)
#save(filtered_correlations, file = 'corr_of_corr_with_SGA/filtered_correlations.RData')

### in the next step of filtering add filtering by positive corr FDR/Bonferroni
queries_with_sig_corr <- filtered_correlations %>%   ## 1225
  filter(grepl('GSP1', query_uniq1) & greater_fdr < 0.05) %>% 
  pull(query_uniq2) %>% unique()
queries_with_two_sig_corr <- filtered_correlations %>% ## 765
  filter(grepl('GSP1', query_uniq1) & greater_fdr < 0.05) %>% 
  group_by(query_uniq2) %>% 
  summarize('sig_count' = n()) %>% 
  filter(sig_count > 1) %>% 
  pull(query_uniq2) %>% unique()
queries_with_bonf_sig_corr <- filtered_correlations %>% ### 501
  filter(grepl('GSP1', query_uniq1) & greater_bonferroni < 0.05) %>% 
  pull(query_uniq2) %>% unique()
queries_with_two_sig_bonf <- filtered_correlations %>% ## 300
  filter(grepl('GSP1', query_uniq1) & greater_bonferroni < 0.05) %>% 
  group_by(query_uniq2) %>% 
  summarize('sig_count' = n()) %>% 
  filter(sig_count > 1) %>% 
  pull(query_uniq2) %>% unique()
# filtered_correlations <- filtered_correlations %>% 
#   filter(query_uniq1 %in% queries_with_sig_corr & query_uniq2 %in% queries_with_sig_corr)
# save(filtered_correlations, file = 'corr_of_corr_with_SGA/filtered_v2_correlations.RData')
# 
# filtered_correlations <- filtered_correlations %>% 
#   filter(query_uniq1 %in% queries_with_two_sig_corr & query_uniq2 %in% queries_with_two_sig_corr)
# save(filtered_correlations, file = 'corr_of_corr_with_SGA/filtered_v3_correlations.RData')

filtered_correlations <- filtered_correlations %>% 
  filter(query_uniq1 %in% queries_with_bonf_sig_corr & query_uniq2 %in% queries_with_bonf_sig_corr) %>% 
  filter(grepl('GSP1', query_uniq1) & ! grepl('GSP1', query_uniq2)) %>% 
  arrange(query_uniq2)
filtered_correlations %>% pull(query_uniq2) %>% unique()
#### final filtered_correlations version 5 has 479 queries
save(filtered_correlations, file = 'corr_of_corr_with_SGA/filtered_v5_correlations.RData')


filtered_correlations <- filtered_correlations %>% 
  filter(query_uniq1 %in% queries_with_two_sig_bonf & query_uniq2 %in% queries_with_two_sig_bonf) %>% 
  filter(grepl('GSP1', query_uniq1) & ! grepl('GSP1', query_uniq2))
filtered_correlations %>% pull(query_uniq2) %>% unique()
#### final filtered_correlations version 6 has 278 queries
save(filtered_correlations, file = 'corr_of_corr_with_SGA/filtered_v6_correlations.RData')

### the versions 2, 3, 4, 5 and 6 filtering still keep it to 22 Gsp1 mutants
filtered_correlations %>% filter(grepl('GSP1', query_uniq1) & (! grepl('GSP1', query_uniq2))) %>% pull(query_uniq1) %>% unique()

### for all the downstream analysis use the version 5 or version 6 - 
#### it's a combination of filtered_correlations (only queries with at least 10 sig interactions with arrays that have sig interactions with at least one Ran mutant)
#### and then additional filtering, that the query needs to have at least one siginificant correlation with any of the mutants (significance based on Boferroni corrected p-value)
