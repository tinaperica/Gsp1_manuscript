### This code makes heatmaps based on correlation data (FDR or Pearson heatmaps) for Figure 1 

library(tidyverse)
library(dendextend)
library(factoextra)
library(cowplot)
source('ucsf_colors.R')
load('Data/filtered_v5_correlations.RData')

pca_order_query_genes <- function(data, dist_method) {
  hc.row <- data %>%
    select(mutant, query, value) %>% 
    spread(query, value) %>% 
    select(-mutant) %>% 
    as.matrix() %>% 
    t() %>% 
    get_dist(method = dist_method) %>% 
    hclust(method = 'single')
  data <- data %>% select(mutant, query, value) %>% spread(mutant, value)
  mat <- as.matrix(data[-1])
  rownames(mat) <- data$query
  row_pcoor <-
    cmdscale(dist(mat), eig = T, k = 1)$point %>%
    as_tibble(rownames = 'gene') %>%
    arrange(`V1`) %>% pull(gene)
  hc.row <- rotate(hc.row, row_pcoor)
  queries_ordered <- hc.row$labels[hc.row$order]
  return(queries_ordered)
}
pca_order_mutants <- function(data, dist_method) {
  hc.row <- data %>%
    select(mutant, query, value) %>% 
    spread(mutant, value) %>% 
    select(-query) %>% 
    as.matrix() %>% 
    t() %>% 
    get_dist(method = dist_method) %>% 
    hclust(method = 'single')
  data <- data %>% select(mutant, query, value) %>% spread(query, value)
  mat <- as.matrix(data[-1])
  rownames(mat) <- data$mutant
  row_pcoor <-
    cmdscale(dist(mat), eig = T, k = 1)$point %>%
    as_tibble(rownames = 'gene') %>%
    arrange(`V1`) %>% pull(gene)
  hc.row <- rotate(hc.row, row_pcoor)
  mutants_ordered <- hc.row$labels[hc.row$order]
  return(rev(mutants_ordered))
}

filtered_correlations_mut <- filtered_correlations %>% 
  filter(grepl('GSP1', query_uniq1) & (! grepl('GSP1', query_uniq2))) %>% 
  select('mutant' = query_uniq1, 'query' = query_uniq2, everything()) %>% 
  mutate('mutant' = substring(mutant, first = 8))
ordered_mutants <- filtered_correlations_mut %>% 
  select(mutant, query, 'value' = greater_fdr) %>% 
  pca_order_mutants(., dist_method = 'pearson')
ordered_queries <- filtered_correlations_mut %>% 
  select(mutant, query, 'value' = greater_fdr) %>% 
  pca_order_query_genes(., dist_method = 'pearson')

make_corr_heatmap <- function(data, value_type = 'pearson') {
  if (value_type == 'pearson') {
    data <- data %>% select(mutant, query, 'value' = pearson)
    plot_limits <- c(-0.3, 0.3)
    fill_breaks <- c(-0.3, 0, 0.3)
    fill_values <- c(ucsf_colors$pink1, 'white', ucsf_colors$green1)
    legend_title <- 'Pearson correlation  '
    tile_color <- 'white'
  } else if (value_type == 'fdr') {
    data <- data %>% select(mutant, query, 'value' = fdr)
    plot_limits <- c(0, 0.1)
    fill_breaks <- c(0.05, 0)
    fill_values <- c(ucsf_colors$purple1, 'white')
    legend_title <- 'FDR of Pearson correlation  '
    tile_color <- ucsf_colors$gray3
  } else if (value_type == 'greater_fdr') {
    data <- data %>% select(mutant, query, 'value' = greater_fdr)
    plot_limits <- c(0, 0.1)
    fill_breaks <- c(0.05, 0)
    fill_values <- c(ucsf_colors$purple1, 'white')
    legend_title <- 'FDR of positive Pearson correlation  '
    tile_color <- ucsf_colors$gray3
  } 
  
  data %>% 
    mutate('mutant' = factor(mutant, ordered_mutants), 
           'query' = factor(query, ordered_queries)) %>% 
    arrange(mutant, query) %>% 
    ggplot(aes(x = query, y = mutant, fill = value)) + 
    geom_tile(color = tile_color, size = 0.05) +
    labs(fill = legend_title) +
    scale_fill_gradientn(breaks = fill_breaks, colours = fill_values, limits = plot_limits, na.value = 'white') +
    theme_classic() +
    ylab('Strong Gsp1 point mutants') +
    xlab('Subset of yeast genes from the SGA dataset with at least one siginificant Pearson correlation with the Gsp1 mutants') +
    theme(
      axis.text.y = element_text(size = 3.5),
      axis.text.x = element_blank(),
      axis.ticks = element_line(size = 0.05),
      axis.ticks.length = unit(0.05, 'cm'),
      axis.line = element_line(size = 0.1),
      axis.title = element_text(size = 4),
      legend.position = 'bottom',
      legend.text = element_text(size = 3.5),
      legend.title = element_text(size = 4),
      legend.key.size = unit(0.18, "cm"),
      legend.margin = margin(0, 0, 0, 0)
    )
}
make_corr_heatmap(filtered_correlations_mut, value_type = 'greater_fdr')
ggsave('E-MAP_Figure1/E-MAP_Figure1/positive_fdr_main_figure.pdf', width = 3.2, height = 1.7)

cor_fig <- make_corr_heatmap(filtered_correlations_mut, value_type = 'pearson')
fdr_fig <- make_corr_heatmap(filtered_correlations_mut, value_type = 'fdr')
cor_sup_fig <- plot_grid(cor_fig, fdr_fig, align = 'h')
ggsave('E-MAP_Figure1/E-MAP_Fig1_supplementary_figures/cor_supplementary_figure.pdf', width = 6.4, height = 2)









