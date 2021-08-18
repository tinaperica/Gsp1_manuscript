library(tidyverse)

# load correlation matrix
load('Data/filtered_v6_correlations.RData')

filtered_correlations <- 
  filtered_correlations %>% 
  filter(!grepl('damp', query_uniq2))


load('Data/filtered_correlations_BNF_2.RData')


queries_v6 <- 
    filtered_correlations %>% 
    pull(query_uniq2) %>% 
    unique()

queries_bnf2 <- 
    filtered_correlations_BNF_2 %>% 
    pull(query_uniq2) %>% 
    unique()

queries_lost <- setdiff(queries_v6, queries_bnf2)

filtered_correlations %>% 
    filter(query_uniq2 %in% queries_lost) %>% 
    filter(greater_p_value >= 0.0001) %>% 
    filter(greater_fdr < 0.05) %>% 
    select(query_uniq2, greater_p_value, greater_fdr) %>% 
    group_by(query_uniq2) %>% 
    summarise(n=n()) %>% 
    filter(n < 2)
    

load('Data/filtered_correlations_BNF_1.RData')
filtered_correlations_BNF_1 %>% 
    filter(query_uniq2 %in% queries_lost) %>% 
    filter(greater_p_value >= 0.0001) %>% 
    filter(greater_fdr < 0.05) %>% 
    select(query_uniq2, greater_p_value, greater_fdr) %>% 
    group_by(query_uniq2) %>% 
    summarise(n=n()) %>% 
    filter(n < 2)



load('Data/spitzemapko_correlations_and_bonferroni_fdr_all.RData')

correlations %>% 
    filter(grepl('GSP1', query_uniq1)) %>% 
    filter(query_uniq2 %in% queries_lost) %>% 
    filter(greater_p_value > 0) %>% 
    filter(greater_p_value < 0.0001)





