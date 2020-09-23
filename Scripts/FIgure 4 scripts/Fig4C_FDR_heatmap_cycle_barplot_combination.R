library(tidyverse)
library(factoextra)
library(cowplot)
library(dendextend)
source('ucsf_colors.R')
load('Data/filtered_v5_correlations.RData')
strong_mutants <- filtered_correlations %>% filter(grepl('GSP1', query_uniq1)) %>% pull(query_uniq1) %>% unique()
mutants_ordered <- rev(read_tsv('Figure4_Multispecificity/Plots/4B_order_of_mutants.txt', col_names = F)$X1)

output_directory <- 'Figure4_Multispecificity/E-MAP_MS_and_biophysics_combined/heatmap_and_barplot_output'
gene_set <- read_tsv('Data/gene_sets.txt') %>% 
  select(query, gene_name, gene_set)
all_queries <- read_tsv('Data/all_queries.txt') %>% 
  select(query, 'gene_name' = name) %>% 
  mutate('gene_set' = 'all')
gene_set <- bind_rows(gene_set, all_queries)
GAP_kinetics <- read_tsv('Data/RanGAP_assay/GAP_kinetics_MichaelisMenten_parameters.txt') %>% 
  select(mutant, 'kcat' = mean_kcat, 'Km' = mean_Km, 
         'kcat_Km' = mean_kcat_Km,	kcat_sd, Km_sd, sd,
         mean_log_kcat, mean_log_Km, log_kcat_sd, log_Km_sd, log_kcat_over_Km, log_sd) %>% 
  unique()
GEF_kinetics <- read_tsv('Data/RanGEF_assay/GEF_kinetics_MichaelisMenten_parameters.txt')
WT_GEF <- GEF_kinetics %>% 
  filter(mutant == "WT")
WT_GAP <- GAP_kinetics %>% 
  filter(mutant == "WT")
GEF_kinetics <- GEF_kinetics %>% 
  mutate('rel_GEF_kcat_Km' = kcat_Km/WT_GEF$kcat_Km) %>% 
  mutate('rel_GEF_kcat' = kcat/WT_GEF$kcat) %>% 
  mutate('rel_GEF_Km' = Km/WT_GEF$Km) %>% 
  select(mutant, rel_GEF_kcat_Km, rel_GEF_kcat, rel_GEF_Km)
GAP_kinetics <- GAP_kinetics %>% 
  mutate('rel_GAP_kcat_Km' = kcat_Km/WT_GAP$kcat_Km) %>% 
  mutate('rel_GAP_kcat' = kcat/WT_GAP$kcat) %>% 
  mutate('rel_GAP_Km' = Km/WT_GAP$Km) %>% 
  select(mutant, rel_GAP_kcat_Km, rel_GAP_kcat, rel_GAP_Km)

kinetics <- GEF_kinetics %>%
  inner_join(., GAP_kinetics, by = 'mutant') %>% 
  mutate('GAP/GEF' = rel_GAP_kcat_Km/rel_GEF_kcat_Km) %>% 
  mutate('GAP/GEF' = ifelse(`GAP/GEF` > 30, 30, `GAP/GEF`)) %>% 
  mutate('rel_GEF_Km' = ifelse(rel_GEF_Km > 30, 30, rel_GEF_Km))

corr <- filtered_correlations %>% 
  filter(query_uniq1 %in% strong_mutants) %>% 
  filter(query_uniq2 %in% unique(gene_set$query)) %>% 
  mutate('query_uniq1' = substring(query_uniq1, first = 8)) %>% 
  filter(query_uniq1 %in% intersect(GAP_kinetics$mutant, GEF_kinetics$mutant)) %>% 
  inner_join(., gene_set, by = c('query_uniq2' = 'query')) %>% 
  select('mutant' = query_uniq1, 'query' = query_uniq2, gene_name, gene_set, greater_fdr, pearson) #%>% 
  #mutate('query' = toupper(query))


gene_sets <- gene_set %>% pull(gene_set) %>% unique()

pca_order_query_genes <- function(data, dist_method) {
  hc.row <- data %>%
    spread(query, greater_fdr) %>% 
    select(-mutant) %>% 
    as.matrix() %>% 
    t() %>% 
    get_dist(method = dist_method) %>% 
    hclust(method = 'single')
  data <- data %>% spread(mutant, greater_fdr)
  mat <- as.matrix(data[-1])
  rownames(mat) <- data$query
  row_pcoor <-
    cmdscale(dist(mat), eig = T, k = 1)$point %>%
    as_tibble(rownames = 'mutant') %>%
    arrange(`V1`) %>% pull(mutant)
  hc.row <- rotate(hc.row, row_pcoor)
  queries_ordered <- hc.row$labels[hc.row$order]
  return(queries_ordered)
}
pca_order_query_genes_by_pearson <- function(data, dist_method) {
  hc.row <- data %>%
    spread(query, pearson) %>% 
    select(-mutant) %>% 
    as.matrix() %>% 
    t() %>% 
    get_dist(method = dist_method) %>% 
    hclust(method = 'single')
  data <- data %>% spread(mutant, pearson)
  mat <- as.matrix(data[-1])
  rownames(mat) <- data$query
  row_pcoor <-
    cmdscale(dist(mat), eig = T, k = 1)$point %>%
    as_tibble(rownames = 'mutant') %>%
    arrange(`V1`) %>% pull(mutant)
  hc.row <- rotate(hc.row, row_pcoor)
  queries_ordered <- hc.row$labels[hc.row$order]
  return(queries_ordered)
}

theme_custom = function(gene_set = 'all') {
  theme_classic() %+replace%
    theme(
      axis.text.x = element_text(colour = ifelse(gene_set == 'all', 'white', 'black'))
    )
}
distance_method <- 'euclidean'

for (i in seq_along(gene_sets)) {
  gene_set_to_plot <- gene_sets[i]
  to_plot <- corr %>% 
    filter(gene_set == gene_set_to_plot) %>% 
    select(mutant, query, greater_fdr)
  to_cluster <- corr %>% 
    filter(gene_set == gene_set_to_plot) %>% 
    select(mutant, query, greater_fdr)
  queries_ordered <- pca_order_query_genes(data = to_cluster, dist_method = distance_method)
  plot <- list()
  plot[[1]] <- to_plot %>% 
    mutate("mutant" = factor(mutant, mutants_ordered), 
           "query" = factor(query, queries_ordered)) %>% 
    arrange(query, mutant) %>% 
    ggplot(aes(y = mutant, x = query, fill = greater_fdr)) + 
    geom_tile(color = ucsf_colors$gray3) + 
    labs(fill = 'FDR of positive Pearson correlation') +
    scale_fill_gradientn(colors = c(ucsf_colors$purple1, 'white'), 
                         limits = c(0, 0.1), breaks=c(0.1, 0), na.value = 'white') +
    theme_classic() +
    ylab('Point mutation in Gsp1') +
    xlab(gene_set_to_plot) +
    theme_custom(gene_set = gene_set_to_plot) +
    theme(text = element_text(size = 6),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = 6),
          axis.ticks = element_line(size = 0.05),
          axis.ticks.length = unit(0.05, 'cm'),
          axis.line = element_line(size = 0.1),
          legend.key.size = unit(0.3, "cm"),
          legend.position = 'bottom',
          legend.box.margin = margin(2,0,0,0)) +
          guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
                 size = guide_legend(title.position = "top", title.hjust = 1))
  kinetics_to_plot <- kinetics %>% 
    filter(mutant %in% mutants_ordered) %>% 
    mutate("mutant" = factor(mutant, mutants_ordered)) %>% 
    arrange(mutant)
  plot[[2]] <- kinetics_to_plot %>% 
    ggplot(aes(mutant, log(`GAP/GEF`), fill = log(`GAP/GEF`))) + 
    geom_bar(stat = "identity", color = 'black', width = 0.7, size = 0.2) + 
    scale_fill_gradientn(colors = c(ucsf_colors$orange1, 'white', ucsf_colors$cyan1), limits = c(-2, 3.5), breaks = c(-2, 0, 3.5)) +
    coord_flip() +
    ylab('relative GAP/GEF efficiency') +
    labs(fill = expression("ln( GAP kcat/Km"[(MUT/WT)]*") / GEF kcat/Km"[(MUT/WT)]*')')) +
    theme_light() +
    xlab(element_blank()) +
    theme(text = element_text(size = 6), 
          axis.title.x = element_text(size = 6, margin = margin(5, 0, 0, 0)), 
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.title = element_text(size = 6),
          axis.ticks = element_line(size = 0.05),
          axis.ticks.length = unit(0.05, 'cm'),
          axis.line = element_line(size = 0.1),
          legend.key.size = unit(0.3, "cm"),
          legend.position = 'bottom',
          legend.box.margin = margin(2,0,0,0),
          legend.text = element_text(size = 3)) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
           size = guide_legend(title.position = "top", title.hjust = 1))
  plot[[3]] <- kinetics_to_plot %>% 
    ggplot(aes(mutant, rel_GEF_kcat_Km, fill = rel_GEF_kcat_Km)) + 
    geom_bar(stat = "identity", color = 'black', width = 0.7, size = 0.2) + 
    scale_fill_gradientn(colors = c('white', ucsf_colors$cyan1), limits = c(0, 1), breaks = c(0, 1)) +
    coord_flip() +
    ylab('relative GEF efficiency') +
    labs(fill = expression("GEF kcat/Km"[(MUT/WT)]*')')) +
    theme_light() +
    xlab(element_blank()) +
    theme(text = element_text(size = 6), 
          axis.title.x = element_text(size = 6, margin = margin(5, 0, 0, 0)), 
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.title = element_text(size = 6),
          axis.ticks = element_line(size = 0.05),
          axis.ticks.length = unit(0.05, 'cm'),
          axis.line = element_line(size = 0.1),
          legend.key.size = unit(0.3, "cm"),
          legend.position = 'bottom',
          legend.box.margin = margin(2,0,0,0),
          legend.text = element_text(size = 3)) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
           size = guide_legend(title.position = "top", title.hjust = 1))
  plot[[4]] <- kinetics_to_plot %>% 
    ggplot(aes(mutant, rel_GAP_kcat_Km, fill = rel_GAP_kcat_Km)) + 
    geom_bar(stat = "identity", color = 'black', width = 0.7, size = 0.2) + 
    scale_fill_gradientn(colors = c('white', ucsf_colors$orange1), limits = c(0, 3.5), breaks = c(0, 3.5)) +
    coord_flip() +
    ylab('relative GAP efficiency') +
    labs(fill = expression("GAP kcat/Km"[(MUT/WT)]*')')) +
    theme_light() +
    xlab(element_blank()) +
    theme(text = element_text(size = 6), 
          axis.title.x = element_text(size = 6, margin = margin(5, 0, 0, 0)), 
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.title = element_text(size = 6),
          axis.ticks = element_line(size = 0.05),
          axis.ticks.length = unit(0.05, 'cm'),
          axis.line = element_line(size = 0.1),
          legend.key.size = unit(0.3, "cm"),
          legend.position = 'bottom',
          legend.box.margin = margin(2,0,0,0),
          legend.text = element_text(size = 3)) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
           size = guide_legend(title.position = "top", title.hjust = 1))
  plot[[5]] <- kinetics_to_plot %>% 
    ggplot(aes(mutant, rel_GEF_Km, fill = rel_GEF_Km)) + 
    geom_bar(stat = "identity", color = 'black', width = 0.7, size = 0.2) + 
    scale_fill_gradientn(colors = c(ucsf_colors$cyan1, 'white'), limits = c(0, 35), breaks = c(0, 35)) +
    coord_flip() +
    ylab(expression('relative GEF K'[m])) +
    labs(fill = expression("relative GEF K"[m(MUT/WT)]*')')) +
    theme_light() +
    xlab(element_blank()) +
    theme(text = element_text(size = 6), 
          axis.title.x = element_text(size = 6, margin = margin(5, 0, 0, 0)), 
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.title = element_text(size = 6),
          axis.ticks = element_line(size = 0.05),
          axis.ticks.length = unit(0.05, 'cm'),
          axis.line = element_line(size = 0.1),
          legend.key.size = unit(0.3, "cm"),
          legend.position = 'bottom',
          legend.box.margin = margin(2,0,0,0),
          legend.text = element_text(size = 3)) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
           size = guide_legend(title.position = "top", title.hjust = 1))
  plot[[6]] <- kinetics_to_plot %>% 
    ggplot(aes(mutant, rel_GAP_Km, fill = rel_GAP_Km)) + 
    geom_bar(stat = "identity", color = 'black', width = 0.7, size = 0.2) + 
    scale_fill_gradientn(colors = c(ucsf_colors$orange1, 'white'), limits = c(0, 7.1), breaks = c(0, 7.1)) +
    coord_flip() +
    ylab(expression('relative GAP K'[m])) +
    labs(fill = expression("relative GAP K"[m(MUT/WT)]*')')) +
    theme_light() +
    xlab(element_blank()) +
    theme(text = element_text(size = 6), 
          axis.title.x = element_text(size = 6, margin = margin(5, 0, 0, 0)), 
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.title = element_text(size = 6),
          axis.ticks = element_line(size = 0.05),
          axis.ticks.length = unit(0.05, 'cm'),
          axis.line = element_line(size = 0.1),
          legend.key.size = unit(0.3, "cm"),
          legend.position = 'bottom',
          legend.box.margin = margin(2,0,0,0),
          legend.text = element_text(size = 3)) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
           size = guide_legend(title.position = "top", title.hjust = 1))
  plot[[7]] <- kinetics_to_plot %>% 
    ggplot(aes(mutant, rel_GEF_kcat, fill = rel_GEF_kcat)) + 
    geom_bar(stat = "identity", color = 'black', width = 0.7, size = 0.2) + 
    scale_fill_gradientn(colors = c('white', ucsf_colors$cyan1), limits = c(0, 3), breaks = c(0, 3)) +
    coord_flip() +
    ylab(expression('relative GEF K'[cat])) +
    labs(fill = expression("relative GEF K"[kcat(MUT/WT)]*')')) +
    theme_light() +
    xlab(element_blank()) +
    theme(text = element_text(size = 6), 
          axis.title.x = element_text(size = 6, margin = margin(5, 0, 0, 0)), 
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.title = element_text(size = 6),
          axis.ticks = element_line(size = 0.05),
          axis.ticks.length = unit(0.05, 'cm'),
          axis.line = element_line(size = 0.1),
          legend.key.size = unit(0.3, "cm"),
          legend.position = 'bottom',
          legend.box.margin = margin(2,0,0,0),
          legend.text = element_text(size = 3)) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
           size = guide_legend(title.position = "top", title.hjust = 1))
  plot[[8]] <- kinetics_to_plot %>% 
    ggplot(aes(mutant, rel_GAP_kcat, fill = rel_GAP_kcat)) + 
    geom_bar(stat = "identity", color = 'black', width = 0.7, size = 0.2) + 
    scale_fill_gradientn(colors = c('white', ucsf_colors$orange1), limits = c(0, 1.5), breaks = c(0, 1.5)) +
    coord_flip() +
    ylab(expression('relative GAP K'[cat])) +
    labs(fill = expression("relative GAP K"[cat(MUT/WT)]*')')) +
    theme_light() +
    xlab(element_blank()) +
    theme(text = element_text(size = 6), 
          axis.title.x = element_text(size = 6, margin = margin(5, 0, 0, 0)), 
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.title = element_text(size = 6),
          axis.ticks = element_line(size = 0.05),
          axis.ticks.length = unit(0.05, 'cm'),
          axis.line = element_line(size = 0.1),
          legend.key.size = unit(0.3, "cm"),
          legend.position = 'bottom',
          legend.box.margin = margin(2,0,0,0),
          legend.text = element_text(size = 3)) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
           size = guide_legend(title.position = "top", title.hjust = 1))
    heatmap_width <- to_plot %>% pull(query) %>% unique() %>% length() %>%  `/`(7.5) #%>% round()
    if (heatmap_width < 1.2) {heatmap_width <- 1.2} else if (heatmap_width < 2) {heatmap_width <- 2}
    if (gene_set_to_plot == 'all') {heatmap_width <- 6}
    barplot_fixed_width <- 0.9
    width_ratio <- heatmap_width/barplot_fixed_width
    file_width <- barplot_fixed_width + heatmap_width
    combn_plots <- plot_grid(plot[[1]], plot[[2]], nrow = 1, align = 'h', rel_widths = c(width_ratio, 1))
    pdf(str_c(output_directory, '/GAP_GEF_ratio/', distance_method, '_and_PCA_dist_hclust_', gene_set_to_plot, '.pdf'), height = 2.7, width = file_width)
    print(combn_plots)
    dev.off()
    combn_plots_separate <- plot_grid(plot[[1]], plot[[3]], plot[[4]], nrow = 1, align = 'h', rel_widths = c(width_ratio, 1, 1))
    pdf(str_c(output_directory, '/GAP_GEF_separate/', distance_method, '_and_PCA_dist_hclust_', gene_set_to_plot, '.pdf'), height = 2.7, width = 1.7*file_width)
    print(combn_plots_separate)
    dev.off()
    combn_plots_GEF <- plot_grid(plot[[1]], plot[[3]], nrow = 1, align = 'h', rel_widths = c(width_ratio, 1))
    pdf(str_c(output_directory, '/GEF_kcat_Km/', distance_method, '_and_PCA_dist_hclust_', gene_set_to_plot, '.pdf'), height = 2.7, width = file_width)
    print(combn_plots_GEF)
    dev.off()
    combn_plots_GAP <- plot_grid(plot[[1]], plot[[4]], nrow = 1, align = 'h', rel_widths = c(width_ratio, 1))
    pdf(str_c(output_directory, '/GAP_kcat_Km/', distance_method, '_and_PCA_dist_hclust_', gene_set_to_plot, '.pdf'), height = 2.7, width = file_width)
    print(combn_plots_GAP)
    dev.off()
    combn_plots_GEF_Km <- plot_grid(plot[[1]], plot[[5]], nrow = 1, align = 'h', rel_widths = c(width_ratio, 1))
    pdf(str_c(output_directory, '/GEF_Km/', distance_method, '_and_PCA_dist_hclust_', gene_set_to_plot, '.pdf'), height = 2.7, width = file_width)
    print(combn_plots_GEF_Km)
    dev.off()
    combn_plots_GAP_Km <- plot_grid(plot[[1]], plot[[6]], nrow = 1, align = 'h', rel_widths = c(width_ratio, 1))
    pdf(str_c(output_directory, '/GAP_Km/', distance_method, '_and_PCA_dist_hclust_', gene_set_to_plot, '.pdf'), height = 2.7, width = file_width)
    print(combn_plots_GAP_Km)
    dev.off()
    combn_plots_GEF_kcat <- plot_grid(plot[[1]], plot[[7]], nrow = 1, align = 'h', rel_widths = c(width_ratio, 1))
    pdf(str_c(output_directory, '/GEF_kcat/', distance_method, '_and_PCA_dist_hclust_', gene_set_to_plot, '.pdf'), height = 2.7, width = file_width)
    print(combn_plots_GEF_kcat)
    dev.off()
    combn_plots_GAP_kcat <- plot_grid(plot[[1]], plot[[8]], nrow = 1, align = 'h', rel_widths = c(width_ratio, 1))
    pdf(str_c(output_directory, '/GAP_kcat/', distance_method, '_and_PCA_dist_hclust_', gene_set_to_plot, '.pdf'), height = 2.7, width = file_width)
    print(combn_plots_GAP_kcat)
    dev.off()
    combn_all_plots <- plot_grid(plot[[1]], plot[[3]], plot[[4]], plot[[5]], plot[[6]], plot[[7]], plot[[8]], nrow = 1, align = 'h', rel_widths = c(width_ratio, 1, 1, 1, 1, 1, 1))
    pdf(str_c(output_directory, '/all_parameters/', distance_method, '_and_PCA_dist_hclust_', gene_set_to_plot, '.pdf'), height = 2.9, width = 2.8*file_width)
    print(combn_all_plots)
    dev.off()
}

main_fig_gene_sets <- c("nuclear pore complex", "spindle assembly checkpoint", "tRNA modification")
plot <- list()
for (i in seq_along(main_fig_gene_sets)) {
  gene_set_to_plot <- main_fig_gene_sets[i]
  to_cluster <- corr %>% 
    filter(gene_set == gene_set_to_plot) %>% 
    select(mutant, query, pearson)
  queries_ordered <- pca_order_query_genes_by_pearson(data = to_cluster, dist_method = distance_method)
  to_plot <- corr %>% 
    filter(gene_set == gene_set_to_plot) %>% 
    select(mutant, query, greater_fdr, pearson) %>% 
    mutate("mutant" = factor(mutant, mutants_ordered), 
           "query" = factor(query, queries_ordered),
           'pearson' = ifelse(pearson <= -0.4, -0.4, pearson),
           'pearson' = ifelse(pearson >= 0.4, 0.4, pearson)) %>% 
    arrange(query, mutant)
    if (i == 1) {
    plot[[gene_set_to_plot]][[1]] <- to_plot %>% 
      ggplot(aes(y = mutant, x = query, fill = greater_fdr)) + 
      geom_tile(color = ucsf_colors$gray3) + 
      labs(fill = 'FDR of positive Pearson correlation') +
      scale_fill_gradientn(colors = c(ucsf_colors$purple1, 'white'), 
                           limits = c(0, 0.1), breaks=c(0.1, 0), na.value = 'white') +
      theme_minimal() +
      ylab('Gsp1 point mutant') +
      xlab(gene_set_to_plot) +
      theme_custom(gene_set = gene_set_to_plot) +
      theme(text = element_text(size = 6),
            axis.title.x = element_text(size = 6),
            axis.title.y = element_text(size = 6),
            axis.text.y = element_text(size = 6),
            axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 6),
            axis.ticks = element_line(size = 0.05),
            axis.ticks.length = unit(0.05, 'cm'),
            axis.line = element_line(size = 0.1),
            legend.key.size = unit(0.3, "cm"),
            legend.position = 'none',
            legend.box.margin = margin(2,0,0,0)) +
      guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
             size = guide_legend(title.position = "top", title.hjust = 1))
    plot[[gene_set_to_plot]][[8]] <- to_plot %>% 
      ggplot(aes(y = mutant, x = query, fill = pearson)) + 
      geom_tile(color = ucsf_colors$gray3) + 
      labs(fill = 'Pearson correlation') +
      scale_fill_gradientn(colors = c(ucsf_colors$pink1, 'white', ucsf_colors$green1),
                        limits = c(-0.4, 0.4), breaks = c(-0.4, 0, 0.4)) +
      theme_minimal() +
      ylab('Gsp1 point mutant') +
      xlab(gene_set_to_plot) +
      theme_custom(gene_set = gene_set_to_plot) +
      theme(text = element_text(size = 6),
            axis.title.x = element_text(size = 6),
            axis.title.y = element_text(size = 6),
            axis.text.y = element_text(size = 6),
            axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 6),
            axis.ticks = element_line(size = 0.05),
            axis.ticks.length = unit(0.05, 'cm'),
            axis.line = element_line(size = 0.1),
            legend.key.size = unit(0.3, "cm"),
            legend.position = 'none',
            legend.box.margin = margin(2,0,0,0)) +
      guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
             size = guide_legend(title.position = "top", title.hjust = 1))
  } else {
    plot[[gene_set_to_plot]][[1]] <- to_plot %>% 
      ggplot(aes(y = mutant, x = query, fill = greater_fdr)) + 
      geom_tile(color = ucsf_colors$gray3) + 
      labs(fill = 'FDR of positive Pearson correlation') +
      scale_fill_gradientn(colors = c(ucsf_colors$purple1, 'white'), 
                           limits = c(0, 0.1), breaks=c(0.1, 0), na.value = 'white') +
      theme_minimal() +
      ylab('') +
      xlab(gene_set_to_plot) +
      theme_custom(gene_set = gene_set_to_plot) +
      theme(text = element_text(size = 6),
            axis.title.x = element_text(size = 6),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 6),
            axis.ticks = element_line(size = 0.05),
            axis.ticks.length = unit(0.05, 'cm'),
            axis.line = element_line(size = 0.1),
            legend.key.size = unit(0.3, "cm"),
            legend.position = 'none',
            legend.box.margin = margin(2,0,0,0)) +
      guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
             size = guide_legend(title.position = "top", title.hjust = 1))
    plot[[gene_set_to_plot]][[8]] <- to_plot %>% 
      ggplot(aes(y = mutant, x = query, fill = pearson)) + 
      geom_tile(color = ucsf_colors$gray3) + 
      labs(fill = 'Pearson correlation') +
      scale_fill_gradientn(colors = c(ucsf_colors$pink1, 'white', ucsf_colors$green1),
                           limits = c(-0.4, 0.4), breaks = c(-0.4, 0, 0.4)) +
      theme_minimal() +
      ylab('') +
      xlab(gene_set_to_plot) +
      theme_custom(gene_set = gene_set_to_plot) +
      theme(text = element_text(size = 6),
            axis.title.x = element_text(size = 6),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 6),
            axis.ticks = element_line(size = 0.05),
            axis.ticks.length = unit(0.05, 'cm'),
            axis.line = element_line(size = 0.1),
            legend.key.size = unit(0.3, "cm"),
            legend.position = 'none',
            legend.box.margin = margin(2,0,0,0)) +
      guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
             size = guide_legend(title.position = "top", title.hjust = 1))
  }
  if (i == 3) {
    kinetics_to_plot <- kinetics %>% 
      filter(mutant %in% mutants_ordered) %>% 
      mutate("mutant" = factor(mutant, mutants_ordered)) %>% 
      arrange(mutant)
    plot[[gene_set_to_plot]][[2]] <- kinetics_to_plot %>% 
      ggplot(aes(mutant, log(`GAP/GEF`), fill = log(`GAP/GEF`))) + 
      geom_bar(stat = "identity", color = 'black', width = 0.7, size = 0.2) + 
      scale_fill_gradientn(colors = c(ucsf_colors$orange1, 'white', ucsf_colors$cyan1), limits = c(-2, 3.5), breaks = c(-2, 0, 3.5)) +
      coord_flip() +
      ylab('relative GAP/GEF efficiency') +
      labs(fill = expression("ln( GAP kcat/Km"[(MUT/WT)]*") / GEF kcat/Km"[(MUT/WT)]*')')) +
      theme_classic() +
      xlab(element_blank()) +
      theme(text = element_text(size = 6), 
            axis.title.x = element_text(size = 6, margin = margin(5, 0, 0, 0)), 
            axis.text.y = element_blank(),
            axis.text.x = element_text(size = 6),
            axis.title = element_text(size = 6),
            axis.ticks = element_line(size = 0.05),
            axis.ticks.length = unit(0.05, 'cm'),
            axis.line = element_line(size = 0.1),
            legend.key.size = unit(0.3, "cm"),
            legend.position = 'none',
            legend.box.margin = margin(2,0,0,0),
            legend.text = element_text(size = 3)) +
      guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
             size = guide_legend(title.position = "top", title.hjust = 1))
    plot[[gene_set_to_plot]][[4]] <- kinetics_to_plot %>% 
      ggplot(aes(mutant, rel_GEF_Km)) + 
      geom_bar(stat = "identity", color = 'black', width = 0.7, size = 0.2) + 
      scale_fill_gradientn(colors = c(ucsf_colors$cyan1, 'white'), limits = c(0, 35), breaks = c(0, 35)) +
      coord_flip() +
      ylab(expression('relative GEF K'[m])) +
      labs(fill = expression("relative GEF K"[m(MUT/WT)]*')')) +
      theme_light() +
      xlab(element_blank()) +
      theme(text = element_text(size = 6), 
            axis.title.x = element_text(size = 6, margin = margin(5, 0, 0, 0)), 
            axis.text.y = element_blank(),
            axis.text.x = element_text(size = 6),
            axis.title = element_text(size = 6),
            axis.ticks = element_line(size = 0.05),
            axis.ticks.length = unit(0.05, 'cm'),
            axis.line = element_line(size = 0.1),
            legend.key.size = unit(0.3, "cm"),
            legend.position = 'none',
            legend.box.margin = margin(2,0,0,0),
            legend.text = element_text(size = 3)) +
      guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
             size = guide_legend(title.position = "top", title.hjust = 1))
    plot[[gene_set_to_plot]][[5]] <- kinetics_to_plot %>% 
      ggplot(aes(mutant, rel_GAP_Km)) + 
      geom_bar(stat = "identity", color = 'black', width = 0.7, size = 0.2) + 
      scale_fill_gradientn(colors = c(ucsf_colors$orange1, 'white'), limits = c(0, 7.1), breaks = c(0, 7.1)) +
      coord_flip() +
      ylab(expression('relative GAP K'[m])) +
      theme_light() +
      xlab(element_blank()) +
      theme(text = element_text(size = 6), 
            axis.title.x = element_text(size = 6, margin = margin(5, 0, 0, 0)), 
            axis.text.y = element_blank(),
            axis.text.x = element_text(size = 6),
            axis.title = element_text(size = 6),
            axis.ticks = element_line(size = 0.05),
            axis.ticks.length = unit(0.05, 'cm'),
            axis.line = element_line(size = 0.1)) +
      guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
             size = guide_legend(title.position = "top", title.hjust = 1))
    plot[[gene_set_to_plot]][[6]] <- kinetics_to_plot %>% 
      ggplot(aes(x = mutant, y = rel_GEF_kcat)) + 
      geom_bar(stat = "identity", color = 'black', width = 0.7, size = 0.2) + 
      scale_fill_gradientn(colors = c('white', ucsf_colors$cyan1), limits = c(0, 3), breaks = c(0, 3)) +
      coord_flip() +
      ylab(expression('relative GEF K'[cat])) +
      labs(fill = expression("relative GEF K"[kcat(MUT/WT)]*')')) +
      theme_light() +
      xlab(element_blank()) +
      theme(text = element_text(size = 6), 
            axis.title.x = element_text(size = 6, margin = margin(5, 0, 0, 0)), 
            axis.text.y = element_blank(),
            axis.text.x = element_text(size = 6),
            axis.title = element_text(size = 6),
            axis.ticks = element_line(size = 0.05),
            axis.ticks.length = unit(0.05, 'cm'),
            axis.line = element_line(size = 0.1),
            legend.key.size = unit(0.3, "cm"),
            legend.position = 'bottom',
            legend.box.margin = margin(2,0,0,0),
            legend.text = element_text(size = 3)) +
      guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
             size = guide_legend(title.position = "top", title.hjust = 1))
    plot[[gene_set_to_plot]][[7]] <- kinetics_to_plot %>% 
      ggplot(aes(mutant, rel_GAP_kcat)) + 
      geom_bar(stat = "identity", color = 'black', width = 0.7, size = 0.2) + 
      scale_fill_gradientn(colors = c('white', ucsf_colors$orange1), limits = c(0, 1.5), breaks = c(0, 1.5)) +
      coord_flip() +
      ylab(expression('relative GAP K'[cat])) +
      labs(fill = expression("relative GAP K"[cat(MUT/WT)]*')')) +
      theme_light() +
      xlab(element_blank()) +
      theme(text = element_text(size = 6), 
            axis.title.x = element_text(size = 6, margin = margin(5, 0, 0, 0)), 
            axis.text.y = element_blank(),
            axis.text.x = element_text(size = 6),
            axis.title = element_text(size = 6),
            axis.ticks = element_line(size = 0.05),
            axis.ticks.length = unit(0.05, 'cm'),
            axis.line = element_line(size = 0.1),
            legend.key.size = unit(0.3, "cm"),
            legend.position = 'bottom',
            legend.box.margin = margin(2,0,0,0),
            legend.text = element_text(size = 3)) +
      guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
             size = guide_legend(title.position = "top", title.hjust = 1))
  }
  heatmap_width <- to_plot %>% pull(query) %>% unique() %>% length() %>%  `/`(6.8) #%>% round()
  if (heatmap_width < 1) {heatmap_width <- 1} else if (heatmap_width < 1.7) {heatmap_width <- 1.7}
  plot[[gene_set_to_plot]][[3]] <- heatmap_width
}

barplot_fixed_width <- 0.9
#width_ratio <- heatmap_width/barplot_fixed_width
file_width <- 1*barplot_fixed_width + plot[["nuclear pore complex"]][[3]] + plot[["spindle assembly checkpoint"]][[3]] + plot[["tRNA modification"]][[3]]
# combn_plots <- plot_grid(plot[["mRNA export"]][[1]], plot[["mRNA export"]][[2]], plot[["spindle assembly checkpoint"]][[1]], plot[["spindle assembly checkpoint"]][[2]], plot[["tRNA modification"]][[1]], plot[["tRNA modification"]][[2]],
#                          nrow = 1, align = 'h', rel_widths = c(plot[["mRNA export"]][[3]]/barplot_fixed_width, 1, plot[["spindle assembly checkpoint"]][[3]]/barplot_fixed_width, 1, plot[["tRNA modification"]][[3]]/barplot_fixed_width, 1))
combn_plots <- plot_grid(plot[["nuclear pore complex"]][[1]], plot[["spindle assembly checkpoint"]][[1]], plot[["tRNA modification"]][[1]], plot[["tRNA modification"]][[2]],
                         nrow = 1, align = 'h', rel_widths = c((1.2*plot[["nuclear pore complex"]][[3]])/barplot_fixed_width, plot[["spindle assembly checkpoint"]][[3]]/barplot_fixed_width, plot[["tRNA modification"]][[3]]/barplot_fixed_width, 1))
pdf(str_c(output_directory, '/', distance_method, '_for_main_figure_2.pdf'), height = 2.2, width = file_width)
print(combn_plots)
dev.off()


barplot_fixed_width <- 0.9
file_width <- 4*barplot_fixed_width + plot[["nuclear pore complex"]][[3]] + plot[["spindle assembly checkpoint"]][[3]] + plot[["tRNA modification"]][[3]]
combn_plots <- plot_grid(plot[["nuclear pore complex"]][[8]], plot[["spindle assembly checkpoint"]][[8]], plot[["tRNA modification"]][[8]], plot[["tRNA modification"]][[4]], plot[["tRNA modification"]][[5]], plot[["tRNA modification"]][[6]], plot[["tRNA modification"]][[7]], 
                         nrow = 1, align = 'h', rel_widths = c(plot[["nuclear pore complex"]][[3]]/barplot_fixed_width, plot[["spindle assembly checkpoint"]][[3]]/barplot_fixed_width, plot[["tRNA modification"]][[3]]/barplot_fixed_width, 1, 1, 1, 1))
pdf(str_c(output_directory, '/', distance_method, '_Ext_Data_Fig_8b.pdf'), height = 2.5, width = 0.8*file_width)
print(combn_plots)
dev.off()



main_fig_gene_sets <- c("Gsp1 partner", "nuclear transport of protein and mRNA", 
                        "CVT pathway mitochondria and vacuole", "actin and polarity", 
                        "transcription regulation and mediator complex", "5' mRNA capping", 
                        'histones and chromatin')
plot <- list()
for (i in seq_along(main_fig_gene_sets)) {
  gene_set_to_plot <- main_fig_gene_sets[i]
  to_plot <- corr %>% 
    filter(gene_set == gene_set_to_plot) %>% 
    select(mutant, query, greater_fdr)
  to_cluster <- corr %>% 
    filter(gene_set == gene_set_to_plot) %>% 
    select(mutant, query, greater_fdr)
  queries_ordered <- pca_order_query_genes(data = to_cluster, dist_method = distance_method)
  if ((i %% 2) != 0) {
    plot[[gene_set_to_plot]][[1]] <- to_plot %>% 
      mutate("mutant" = factor(mutant, mutants_ordered), 
             "query" = factor(query, queries_ordered)) %>% 
      arrange(query, mutant) %>% 
      ggplot(aes(y = mutant, x = query, fill = greater_fdr)) + 
      geom_tile(color = ucsf_colors$gray3) + 
      labs(fill = 'FDR of positive Pearson correlation') +
      scale_fill_gradientn(colors = c(ucsf_colors$purple1, 'white'), 
                           limits = c(0, 0.1), breaks=c(0.1, 0), na.value = 'white') +
      theme_minimal() +
      ylab('Gsp1 point mutant') +
      xlab(gene_set_to_plot) +
      theme_custom(gene_set = gene_set_to_plot) +
      theme(text = element_text(size = 6),
            axis.title.x = element_text(size = 6),
            axis.title.y = element_text(size = 6),
            axis.text.y = element_text(size = 6),
            axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 6),
            axis.ticks = element_line(size = 0.05),
            axis.ticks.length = unit(0.05, 'cm'),
            axis.line = element_line(size = 0.1),
            legend.key.size = unit(0.3, "cm"),
            legend.position = 'bottom',
            legend.box.margin = margin(2,0,0,0)) +
      guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
             size = guide_legend(title.position = "top", title.hjust = 1))
  } else {
    plot[[gene_set_to_plot]][[1]] <- to_plot %>% 
      mutate("mutant" = factor(mutant, mutants_ordered), 
             "query" = factor(query, queries_ordered)) %>% 
      arrange(query, mutant) %>% 
      ggplot(aes(y = mutant, x = query, fill = greater_fdr)) + 
      geom_tile(color = ucsf_colors$gray3) + 
      labs(fill = 'FDR of positive Pearson correlation') +
      scale_fill_gradientn(colors = c(ucsf_colors$purple1, 'white'), 
                           limits = c(0, 0.1), breaks=c(0.1, 0), na.value = 'white') +
      theme_minimal() +
      ylab('') +
      xlab(gene_set_to_plot) +
      theme_custom(gene_set = gene_set_to_plot) +
      theme(text = element_text(size = 6),
            axis.title.x = element_text(size = 6),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = 6),
            axis.ticks = element_line(size = 0.05),
            axis.ticks.length = unit(0.05, 'cm'),
            axis.line = element_line(size = 0.1),
            legend.key.size = unit(0.3, "cm"),
            legend.position = 'bottom',
            legend.box.margin = margin(2,0,0,0)) +
      guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
             size = guide_legend(title.position = "top", title.hjust = 1))
  }
  if ((i %% 2) == 0 | i == length(main_fig_gene_sets)) {
    kinetics_to_plot <- kinetics %>% 
      filter(mutant %in% mutants_ordered) %>% 
      mutate("mutant" = factor(mutant, mutants_ordered)) %>% 
      arrange(mutant)
    plot[[gene_set_to_plot]][[2]] <- kinetics_to_plot %>% 
      ggplot(aes(mutant, log(`GAP/GEF`), fill = log(`GAP/GEF`))) + 
      geom_bar(stat = "identity", color = 'black', width = 0.7, size = 0.2) + 
      scale_fill_gradientn(colors = c(ucsf_colors$orange1, 'white', ucsf_colors$cyan1), limits = c(-2, 3.5), breaks = c(-2, 0, 3.5)) +
      coord_flip() +
      ylab('relative GAP/GEF efficiency') +
      labs(fill = expression("ln( GAP kcat/Km"[(MUT/WT)]*") / GEF kcat/Km"[(MUT/WT)]*')')) +
      theme_classic() +
      xlab(element_blank()) +
      theme(text = element_text(size = 6), 
            axis.title.x = element_text(size = 6, margin = margin(5, 0, 0, 0)), 
            axis.text.y = element_blank(),
            axis.text.x = element_text(size = 6),
            axis.title = element_text(size = 6),
            axis.ticks = element_line(size = 0.05),
            axis.ticks.length = unit(0.05, 'cm'),
            axis.line = element_line(size = 0.1),
            legend.key.size = unit(0.3, "cm"),
            legend.position = 'bottom',
            legend.box.margin = margin(2,0,0,0),
            legend.text = element_text(size = 3)) +
      guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
             size = guide_legend(title.position = "top", title.hjust = 1))
  }
  heatmap_width <- to_plot %>% pull(query) %>% unique() %>% length() %>%  `/`(7.5) #%>% round()
  if (heatmap_width < 1.2) {heatmap_width <- 1.2} else if (heatmap_width < 2) {heatmap_width <- 2}
  plot[[gene_set_to_plot]][[3]] <- heatmap_width
}

barplot_fixed_width <- 0.9
width_ratio <- heatmap_width/barplot_fixed_width
file_width <- 1*barplot_fixed_width + plot[["Gsp1 partner"]][[3]] + plot[["nuclear transport of protein and mRNA"]][[3]]
combn_plots <- plot_grid(plot[["Gsp1 partner"]][[1]], plot[["nuclear transport of protein and mRNA"]][[1]], plot[["nuclear transport of protein and mRNA"]][[2]],
                         nrow = 1, align = 'h', rel_widths = c(plot[["Gsp1 partner"]][[3]]/barplot_fixed_width, plot[["nuclear transport of protein and mRNA"]][[3]]/barplot_fixed_width, 1))
pdf(str_c(output_directory, '/', distance_method, '_Ext_Fig9a.pdf'), height = 2.9, width = file_width*1.1)
print(combn_plots)
dev.off()


file_width <- 1*barplot_fixed_width + plot[["transcription regulation and mediator complex"]][[3]] + plot[["5' mRNA capping"]][[3]]
combn_plots <- plot_grid(plot[["transcription regulation and mediator complex"]][[1]], plot[["5' mRNA capping"]][[1]], plot[["5' mRNA capping"]][[2]],
                         nrow = 1, align = 'h', rel_widths = c(plot[["transcription regulation and mediator complex"]][[3]]/barplot_fixed_width, plot[["5' mRNA capping"]][[3]]/barplot_fixed_width, 1))
pdf(str_c(output_directory, '/', distance_method, '_Ext_Fig9b.pdf'), height = 2.9, width = file_width*1.1)
print(combn_plots)
dev.off()

file_width <- 1*barplot_fixed_width + plot[["CVT pathway mitochondria and vacuole"]][[3]] + plot[["actin and polarity"]][[3]]
combn_plots <- plot_grid(plot[["CVT pathway mitochondria and vacuole"]][[1]], plot[["actin and polarity"]][[1]], plot[["actin and polarity"]][[2]],
                         nrow = 1, align = 'h', rel_widths = c(plot[["CVT pathway mitochondria and vacuole"]][[3]]/barplot_fixed_width, plot[["actin and polarity"]][[3]]/barplot_fixed_width, 1))
pdf(str_c(output_directory, '/', distance_method, '_Ext_Fig9c.pdf'), height = 2.9, width = file_width*1.1)
print(combn_plots)
dev.off()


file_width <- 1*barplot_fixed_width + plot[["histones and chromatin"]][[3]]
combn_plots <- plot_grid(plot[["histones and chromatin"]][[1]], plot[["histones and chromatin"]][[2]],
                         nrow = 1, align = 'h', rel_widths = c(plot[["histones and chromatin"]][[3]]/barplot_fixed_width, 1))
pdf(str_c(output_directory, '/', distance_method, '_Ext_Fig9d.pdf'), height = 2.9, width = file_width*1.1)
print(combn_plots)
dev.off()

