library(tidyverse)
library(cowplot)

# load correlation matrix
load('Data/filtered_v6_correlations.RData')
filtered_correlations <- filter(filtered_correlations, !query_uniq2 %in% c('crm1_damp','yrb2_damp'))

# clustering function for heatmaps
clustfn <- function(mat, dist_method, clust_method) {
  
  # use euclidean clustering to get first principal coordinate ordering
  pcoor <-
    cmdscale(get_dist(mat, method = 'euclidean'), eig = T, k = 1)$point %>%
    as_tibble(rownames = 'var') %>%
    arrange(`V1`) %>%
    pull(var)
  
  # use the distance methods and clustering method of choice to cluster
  hc <-
    mat %>%
    get_dist(method = dist_method) %>%
    hclust(method = clust_method) %>%
    rotate(pcoor)
  return(hc)
}

group1 <- c('D79S','T34Q','T34E','K101R','D79A','T34G')
group2 <- c('T34A','Q147E','R108I','R108L','G80A','Y157A','H141E')
group3 <- c('H141R','R108Y','R108Q','R108G','Y148I','H141I','R112A','R112S','R78K')

set.seed(3)
holdout1 <- sample(group1,1)
holdout2 <- sample(group2,1)
holdout3 <- sample(group3,1)

df_holdouts <- filtered_correlations %>%
  filter(grepl('GSP1', query_uniq1) & (! grepl('GSP1', query_uniq2))) %>%
  mutate('mutant' = substring(query_uniq1, first = 8)) %>%
  select(mutant, 'query' = query_uniq2, greater_fdr) %>% 
  filter(mutant %in% c(holdout1, holdout2, holdout3)) %>% 
  rename('group'=mutant)

# transform the correlations table into a wide format matrix, rows are mutants, columns are queries
df_corr <-
  filtered_correlations %>%
  filter(grepl('GSP1', query_uniq1) & (! grepl('GSP1', query_uniq2))) %>%
  mutate('mutant' = substring(query_uniq1, first = 8)) %>%
  select(mutant, 'query' = query_uniq2, greater_fdr) %>% 
  mutate(group = case_when(
             mutant %in% group1 ~ 'I',
             mutant %in% group2 ~ 'II',
             mutant %in% group3 ~ 'III')) %>% 
  group_by(group, query) %>% 
  summarise(greater_fdr = mean(greater_fdr)) %>% 
  bind_rows(df_holdouts) %>% 
  spread(group, greater_fdr)
  
theme <- theme_light() +
    theme(
        text = element_text(family = "Helvetica", size = 6),
        axis.title = element_text(size = 6), axis.text = element_text(size = 6),
        axis.ticks = element_line(size = 0.05), axis.ticks.length = unit(0.05, 'cm'),
        strip.text.x = element_text(size = 6)
    )


p11 <- ggplot(df_corr, aes_string(x=holdout1,y='I')) + geom_point() + theme
p12 <- ggplot(df_corr, aes_string(x=holdout1,y='II')) + geom_point() + theme
p13 <- ggplot(df_corr, aes_string(x=holdout1,y='III')) + geom_point() + theme

p21 <- ggplot(df_corr, aes_string(x=holdout2,y='I')) + geom_point() + theme
p22 <- ggplot(df_corr, aes_string(x=holdout2,y='II')) + geom_point() + theme
p23 <- ggplot(df_corr, aes_string(x=holdout2,y='III')) + geom_point() + theme

p31 <- ggplot(df_corr, aes_string(x=holdout3,y='I')) + geom_point() + theme
p32 <- ggplot(df_corr, aes_string(x=holdout3,y='II')) + geom_point() + theme
p33 <- ggplot(df_corr, aes_string(x=holdout3,y='III')) + geom_point() + theme

p <- cowplot::plot_grid(p11,p12,p13,p21,p22,p23,p31,p32,p33,
        nrow=3
)

title <- ggdraw() + 
  draw_label(
    "Pairwise comparisons of holdouts from each group\nand the mean p-value of the group",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

plot_grid(title, p, ncol=1, rel_heights = c(0.1, 1))




## Can we do this systematically for all 

group1 <- c('D79S','T34Q','T34E','K101R','D79A','T34G')
group2 <- c('T34A','Q147E','R108I','R108L','G80A','Y157A','H141E')
group3 <- c('H141R','R108Y','R108Q','R108G','Y148I','H141I','R112A','R112S','R78K')
mutants <- c(group1, group2, group3)

df <-
  filtered_correlations %>%
  filter(grepl('GSP1', query_uniq1) & (! grepl('GSP1', query_uniq2))) %>%
  mutate('mutant' = substring(query_uniq1, first = 8)) %>%
  select(mutant, 'query' = query_uniq2, greater_fdr)


holdout_corrs <- matrix(nrow=length(mutants), ncol=3)
rownames(holdout_corrs) <- mutants
colnames(holdout_corrs) <- c('I','II','III')

for (i in seq(1,length(mutants))) {
    m <- mutants[i]
    df_mutant <- df %>% 
        filter(mutant==m) %>% 
        mutate(group='holdout') %>% 
        select(-mutant)
    df_tmp <- df %>% 
        filter(mutant!=m) %>% 
        mutate(group = case_when(
             mutant %in% group1 ~ 'I',
             mutant %in% group2 ~ 'II',
             mutant %in% group3 ~ 'III')) %>% 
        group_by(group, query) %>% 
        summarise(greater_fdr = mean(greater_fdr), .groups='drop') %>% 
        bind_rows(df_mutant) %>% 
        spread(group, greater_fdr)

    holdout_corrs[m,'I'] <- cor(df_tmp$holdout, df_tmp$I, method='pearson')
    holdout_corrs[m,'II'] <- cor(df_tmp$holdout, df_tmp$II, method='pearson')
    holdout_corrs[m,'III'] <- cor(df_tmp$holdout, df_tmp$III, method='pearson')

    # Ward's criterion (squared euclidean distance)    
    # holdout_corrs[m,'I'] <- sum((df_tmp$holdout - df_tmp$I) ^ 2)
    # holdout_corrs[m,'II'] <- sum((df_tmp$holdout - df_tmp$II) ^ 2)
    # holdout_corrs[m,'III'] <- sum((df_tmp$holdout - df_tmp$III) ^ 2)
    
}



plot_corrs <- function(og) {
  holdout_corrs %>% 
  as_tibble(rownames='mutant') %>% 
  pivot_longer(cols=c('I','II','III'), names_to='group',values_to='cor') %>% 
  mutate(original_group = case_when(
             mutant %in% group1 ~ 'I',
             mutant %in% group2 ~ 'II',
             mutant %in% group3 ~ 'III')) %>% 
  mutate(alpha = ifelse(group==og, 1, 0)) %>% 
  mutate(mutant = factor(mutant, levels=mutants)) %>% 
  filter(original_group == og) %>% 
  ggplot(aes(x=group, y=cor, alpha=alpha, group=original_group)) +
  geom_bar(position='dodge', stat='identity', color='black', size=0.2) +
  facet_wrap(~mutant, nrow=3) +
  ylim(c(1.1*min(holdout_corrs), 1)) +
  theme +
  theme(legend.position='none')
}

p1 <- plot_corrs('I') + xlab('') + ylab('Pearson correlation with centroid of cluster')
p2 <- plot_corrs('II') + xlab('Mutant group') + ylab('')
p3 <- plot_corrs('III') + xlab('') + ylab('')
plot_grid(p1,p2,p3, nrow=1, rel_widths = c(2,3,3))
ggsave('Revisions2/Holdouts_for_reviewer_2.pdf', width = 4.5, height = 2)

holdout_corrs['G80A',]
