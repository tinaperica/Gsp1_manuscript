library(tidyverse)
library(gmp) #### for factorize

### these are bad library strains in the Krogan chromatin E-MAP library
bad_strains <- c('YBR084C-A', 'YBR098W', 'YCR036W',  'YCR044C',  'YCR077C',  'YHR090C', 
                 'YJL117W', 'YJL190C', 'YKR024C', 'YKR062W', 'YML008C', 'YML051W', 'YMR231W', 
                 'YNL330C', 'YPR072W', 'YDL192W', 'YDL082W', 'YDR443C', 'YDL135C', 'YDL135C')
load("Data/spitzemapko.rda")

### spitzemapko has the damp queries removed, and it also keeps only ko array genes 
### (so we can compare them to our ko genes)
### damp queries are removed because the Costanzo paper says 
### that for many of the damps they couldn't confirm lower expression
### however, keep damps for CRM1 and YRB2, because those are the only queries we have for those Gsp1 partners
#### partners for which there is no data in the Costanzo dataset: NMD5, NUP116, KAP104, KAP120, (PSE1 and MTR10 are only screened against the TSA library) 
#load("Data/spitzemap.rda")
load("Data/SGA_Scaling/SGA_2016_full_scaled.rda")
damps_to_add <- SGA_scaled_to_EMAP %>% 
  filter(query_allele_name %in% c('crm1_damp', 'yrb2_damp')) %>% 
  filter(arraytype_temp == 'DMA30') %>% #### using this array type because it is much bigger: 3349 values vs 890
  select(names(spitzemapko)) %>% 
  arrange(query_allele_name, array_allele_name)

#### distribution of scores
damps_to_add %>% ggplot(aes(x = score, color = query_allele_name)) + geom_density()
## looks like all of them have significant score values, might mean they are real damps

#### how are their correlations with our Gsp1 mutants (we should see positive correlations as they are partners)
spitzemapko %>% filter(query_allele_name %in% c('GSP1 - T34E', 'GSP1 - R108L')) %>% 
  inner_join(damps_to_add, by = 'array_ORF') %>% 
  ggplot(aes(x = score.x, y = score.y)) + geom_point() + facet_grid(query_allele_name.x ~ query_allele_name.y)
### they all seem to have a phenotype so I'll keep them
spitzemapko <- bind_rows(spitzemapko, damps_to_add)
sga_index <- spitzemapko %>% select(query_allele_name, query_ORF) %>% unique()
write_tsv(sga_index, 'Data/E-MAP/spitzemapko_query_orf_index.txt')
array_ORFs_from_SGA <- spitzemapko %>% 
  filter(interaction_network != "gsp1_pEMAP") %>% 
  pull(array_ORF) %>% unique()
array_orfs_to_keep <- spitzemapko %>% 
  filter(interaction_network == "gsp1_pEMAP" & ! (array_ORF %in% bad_strains) & array_ORF %in% array_ORFs_from_SGA) %>% 
  pull(array_ORF) %>% unique()

spitzemapko_for_corr <- spitzemapko %>% 
  filter(array_ORF %in% array_orfs_to_keep) %>% 
  select(query_allele_name, query_ORF, array_ORF, score)

all_queries <- spitzemapko_for_corr %>% pull(query_allele_name) %>% unique()

###### remove query_strain_ids with too many NA
queries_to_keep <- spitzemapko_for_corr %>% 
  filter(! is.na(score)) %>% 
  group_by(query_allele_name) %>% 
  summarise("count_array" = n()) %>% 
  arrange(count_array) 
queries_to_keep %>% filter(grepl('rna1', query_allele_name) | 
                             grepl('yrb', query_allele_name) |
                           grepl('srm1', query_allele_name) | 
                             grepl('ntf', query_allele_name) |
                           grepl('kap', query_allele_name) |
                           grepl('crm1', query_allele_name))

queries_to_keep %>% 
  ggplot(aes(x = count_array)) + geom_density()
queries_to_keep <- queries_to_keep %>% 
  filter(count_array > 950)
spitzemapko_for_corr <- spitzemapko_for_corr %>% 
  filter(query_allele_name %in% queries_to_keep$query_allele_name)

#### How many (non-Gsp1) query alleles are there in spitzemapko_for_correlations?
spitzemapko_for_corr %>% filter(! query_ORF == 'YLR293C') %>% pull(query_allele_name) %>% unique() %>% length()
### 3608
### How many unique ORFs is that?
spitzemapko_for_corr %>% filter(! query_ORF == 'YLR293C') %>% pull(query_ORF) %>% unique() %>% length()
### 3369

### add a weight to the score based on a score confidence spline
# The fit spline produces the following values:
#   -6    -5    -4    -3    -2    -1    0     1     2     3     4     5     6
# 1.00  0.99  0.95  0.83  0.43  0.14  0.02  0.07  0.26  0.65  0.86  0.95  0.99
# 
# Note that the spline fits S-score < -5.6 to values slightly (by 1e-3) above 1.
# If these confidence values are being used as weights for scores, simply cap all
# weights at 1.0.
d <- data.frame(avg_S_score = seq(-6, 6, by = 0.5),
                confidence = c(1, 1, 0.99, 0.98, 0.95, 0.9, 0.85, 0.65, 0.45,
                               0.2, 0.15, 0.08, 0, 0.05, 0.07, 0.15,
                               0.22, 0.5, 0.65, 0.75, 0.87, 0.93, .95, .98, 0.99))
confidence_spline <- smooth.spline(x = d$avg_S_score, y = d$confidence, spar = 0.3)
predict(confidence_spline, -2)
spitzemapko_for_corr <- spitzemapko_for_corr %>% 
  mutate("weight" = predict(confidence_spline, score)$y) %>% 
  mutate("weight" = ifelse(weight > 1, 1, weight))


### now make task.info files that will define which pairwise correlations are calculated in 
## that particular task on the cluster
# a general function that makes task files when pairs and steps are defined
make_task_files <- function(pairs, step, task_outpath) {
  tasks <- seq(1, length(pairs)/2, step)
  for (t in seq_along(tasks)) {
    first_pair <- tasks[t]
    last_pair <- first_pair + step - 1
    task.info <- list()
    task.info[["pairs"]] <- pairs[, first_pair:last_pair]
    outfilename <- str_c(task_outpath, first_pair, "_task_info.RData", sep = "")
    save(task.info, file = outfilename)
  }
}

mutants <- spitzemapko_for_corr %>% 
  filter(grepl(query_allele_name, pattern = "GSP1 - ")) %>%
  filter(! query_allele_name == "GSP1 - NTER3XFLAG WT") %>% 
  filter(! query_allele_name == "GSP1 - CTER3XFLAG WT") %>% 
  filter(! query_allele_name == "GSP1 - T34N") %>% 
  pull(query_allele_name) %>% unique()

all_genes_and_mutants <- spitzemapko_for_corr %>%
  filter(! query_allele_name == "GSP1 - NTER3XFLAG WT") %>% 
  filter(! query_allele_name == "GSP1 - CTER3XFLAG WT") %>% 
  filter(! query_allele_name == "GSP1 - T34N") %>% 
  pull(query_allele_name) %>% unique()
  
length(all_genes_and_mutants)
genes <- all_genes_and_mutants[! all_genes_and_mutants %in% mutants]
length(genes)

all_gsp1_pairs <- expand.grid('mutants' = mutants, 'genes' = genes)
all_gsp1_pairs <- tibble('mutants' = as.character(all_gsp1_pairs$mutants), 'genes' = as.character(all_gsp1_pairs$genes)) 
### all the pairs, including the mutants, to calculate all correlations/similarities
#all_pairs <- combn(all_genes_and_mutants, 2) ### all pairs (including mutants)
#pairs_to_save <- tibble("gene1" = all_pairs[1, ], "gene2" = all_pairs[2, ]) 
write_tsv(all_gsp1_pairs, "Data/pairs.txt")

### make task files for all pairs
(n_pairs <- length(all_pairs)/2)
(factors <- factorize(n_pairs))
(step <- 47 * 13 * 3)
n_pairs / step
### 1-6721611:1833 for correlations - 3667 jobs
#make_task_files(all_pairs, step, "Data/all_correlations_task_info/")

spitzemapko_for_corr <- spitzemapko_for_corr %>% 
  select(-query_ORF)
#### save spitzemapko for correlations
save(spitzemapko_for_corr, file = "Data/spitzemapko_for_corr.rda")
