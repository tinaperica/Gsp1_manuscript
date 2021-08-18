library(tidyverse)
library(factoextra)
library(cowplot)
library(dendextend)
library(scales)
library(ggnewscale)
set.seed(2017)
source('ucsf_colors.R')
### output directory
output_directory <- 'Revisions2/Fig_and_analysis_for_StatReviewer/nr_corr/heatmap_and_barplot_output'
### input files
load('Data/filtered_v5_correlations.RData')

filtered_correlations <- 
  filtered_correlations %>% 
  filter(!grepl('damp', query_uniq2)
  
strong_mutants <- filtered_correlations %>% filter(grepl('GSP1', query_uniq1)) %>% pull(query_uniq1) %>% unique()
mutants_ordered <- read_tsv('Data/4B_order_of_mutants.txt', col_names = F)$X1
gene_set_data <- read_tsv('Data/gene_sets.txt') %>% 
  select(query, gene_name, gene_set)

all_queries <- read_tsv('Data/all_queries.txt') %>% 
  select(query, 'gene_name' = name) %>% 
  mutate('gene_set' = 'all')
gene_set_data <- bind_rows(gene_set_data, all_queries)
gene_sets <- gene_set_data %>% pull(gene_set) %>% unique()

# This is no longer capped
kinetics <- read_tsv("Data/kinetics_data_relative_to_WT.txt") %>% 
  filter(! measure %in% c('int')) %>% 
  filter(! grepl('gamma', measure))

# only plot mutants that have both GAP and GEF kinetic measurements
mutants_with_complete_kinetics <-
  kinetics %>%
  select(mutant, measure, rel_to_WT) %>% 
  pivot_wider(values_from = rel_to_WT, names_from = 'measure') %>% 
  drop_na(GAP_kcat_Km, GEF_kcat_Km) %>% 
  pull(mutant) %>%
  unique()

# This is the ordering for plots
mutants_for_figures <- intersect(mutants_ordered, mutants_with_complete_kinetics)

#### This is the part where I make the correlations non-redundant
#### nr_corr data has a randomly selected allele for those genes that have multiple alleles
### e.g. multiple different temperature sensitive mutants, in the spitzemap
nr_alleles <- filtered_correlations %>% 
  filter(query_uniq2 %in% unique(gene_set$query)) %>% 
  inner_join(., gene_set, by = c('query_uniq2' = 'query')) %>% 
  select(query_uniq2, gene_name) %>% 
  unique() %>% 
  group_by(gene_name) %>% 
  sample_n(1) %>% 
  ungroup() %>% 
  pull(query_uniq2)
corr <- filtered_correlations %>% 
  filter(query_uniq1 %in% strong_mutants) %>% 
  filter(query_uniq2 %in% unique(gene_set_data$query)) %>% 
  mutate('query_uniq1' = substring(query_uniq1, first = 8)) %>% 
  filter(query_uniq1 %in% mutants_for_figures) %>% 
  inner_join(., gene_set_data, by = c('query_uniq2' = 'query')) %>% 
  select('mutant' = query_uniq1, 'query' = query_uniq2, gene_name, gene_set, greater_fdr, pearson) %>% 
  filter(query %in% nr_alleles)



### PLOT KINETICS
# mod_plot <- kinetics %>% 
#   filter(mutant %in% mutants_for_figures) %>% 
#   mutate('mutant' = factor(mutant, mutants_for_figures)) %>% 
#   filter(measure %in% c('GAP_kcat_Km', 'GEF_kcat_Km')) %>% 
#   ggplot(aes(x = mutant, y = ln_rel_to_WT, color = measure)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_point(size = 3, alpha = 0.6, shape = 15) +
#   scale_color_manual(values = c(ucsf_colors$orange1, ucsf_colors$cyan1)) +
#   theme_light() +
#   scale_y_continuous(position = 'left', limits = c(-3, 1.2)) +
#   xlab('\nGsp1 point mutant') +
#   ylab('Relative GAP or GEF efficiency\n') +
#   theme(text = element_text(size = 6),
#         axis.title.x = element_text(size = 6, margin = margin(0, 0, 0, 0)), 
#         axis.text.y = element_text(size = 6),
#         axis.text.x = element_text(size = 6, angle = 90),
#         axis.title = element_text(size = 6),
#         axis.ticks = element_line(size = 0.05),
#         axis.ticks.length = unit(0.05, 'cm'),
#         axis.line = element_line(size = 0.1),
#         legend.key.size = unit(0.3, "cm"),
#         legend.position = 'right',
#         legend.box.margin = margin(0,0,0,0),
#         legend.text = element_text(size = 3))
# 

### START LINE PLOT TRY 2020-06-30
mut_groups <- list('I' = c('D79S', 'T34Q', 'T34E', 'K101R', 'T34G', 'D79A'),
                   'II' = c('T34A', 'Q147E', 'R108I', 'R108L', 'G80A', 'Y157A', 'H141E'),
                   'III' = c('H141R', 'R108Y', 'R108Q', 'R108G', 'Y148I', 'H141I', 'R112A', 'R112S', 'R78K')) 
outliers <- c('K101R', 'R112S', 'R78K')
lp_df <-
  kinetics %>% 
  filter(mutant %in% mutants_for_figures) %>% 
  mutate('mutant' = factor(mutant, mutants_for_figures)) %>% 
  filter(measure %in% c('GAP_kcat_Km', 'GEF_kcat_Km')) %>% 
  mutate('group' = case_when(mutant %in% mut_groups$I ~ 'I',
                             mutant %in% mut_groups$II ~ 'II',
                             mutant %in% mut_groups$III ~ 'III')) %>%
  mutate(ln_rel_to_WT = ifelse(ln_rel_to_WT < -3, -3, ln_rel_to_WT)) %>% 
  mutate(linetype = ifelse(mutant %in% outliers, 'outlier', 'non-outlier')) %>% 
  mutate(shape = ifelse(mutant %in% outliers, 'outlier', 'non-outlier'))
line_plot <- ggplot(data = lp_df,
       mapping = aes(x = measure, y = ln_rel_to_WT, color = measure, group = mutant)) +
  facet_grid(cols = vars(group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = 'gray') +
  geom_point(aes(shape = shape), size = 1.5, alpha = 0.6, stroke = .5) +
  geom_line(aes(linetype = linetype), color = ucsf_colors$gray1) +
  stat_summary(data = filter(lp_df, mutant != 'K101R'),
               mapping=aes(group=measure), geom = "crossbar",
               fun.y = median, fun.ymin = median, fun.ymax = median,
               fatten = 0, width = 0.5, color = ucsf_colors$pink1) + 
  scale_color_manual(values = c(ucsf_colors$orange1, ucsf_colors$cyan1)) +
  scale_linetype_manual(limits = c('non-outlier', 'outlier'),
                        values = c('solid', 'dashed'),
                        labels = c('non-outlier', 'outliers: K101R, R112S, R78K'),
  ) +
  scale_shape_manual(limits = c('non-outlier', 'outlier'),
                     values = c(19, 1),
                     labels = c('non-outlier', 'outliers: K101R, R112S, R78K')) +
  theme_light() +
  # scale_y_continuous(position = 'left', limits = c(-3, 1.2)) +
  xlab('\nGsp1 point mutant') +
  ylab('Relative GAP or GEF efficiency\n') +
  theme(text = element_text(size = 6),
        axis.title.x = element_text(size = 6, margin = margin(0, 0, 0, 0)), 
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6, angle = 90),
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = 0.05),
        axis.ticks.length = unit(0.05, 'cm'),
        axis.line = element_line(size = 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.position = 'right',
        legend.box.margin = margin(0,0,0,0),
        legend.text = element_text(size = 3))
#ggsave('Revisions/Main Figures/Figure4/Fig4_LinkedPoints_LOG.pdf', height = 5, width = 5)
#dev.off()
### END LINE PLOT TRY 2020-06-30


### SINA PLOT FOR Fig. 4B - print it together with Fig. 4c so they are nicely aligned
# clean correlations dataset, so each row is a correlation between a mutant and a strain
mutant_group_order = c('I', 'II', 'III')
mutant_groups_for_sina <- 
  list('I' = c('D79S', 'T34Q', 'T34E', 'K101R', 'T34G', 'D79A'),
       'II' = c('T34A', 'Q147E', 'R108I', 'R108L', 'G80A', 'Y157A', 'H141E'),
       'III' = c('H141R', 'R108Y', 'R108Q', 'R108G', 'Y148I', 'H141I', 'R112A', 'R112S', 'R78K')) %>% 
  stack() %>% 
  rename('mutant' = 'values', 'mutant_group' = 'ind') %>% 
  mutate(mutant_group = factor(mutant_group, levels = mutant_group_order))

corr_for_sina <- corr %>%
  filter(mutant %in% unlist(mutant_groups_for_sina)) %>% 
  select(mutant, 'strain' = query, pearson, greater_fdr) %>%
  left_join(mutant_groups_for_sina, by = 'mutant') 
# annotate strains based on gene set  
gene_sets_order <- c('nuclear pore complex', 'spindle assembly checkpoint', 'tRNA modification', 'other')
gene_sets_colors <- c(ucsf_colors$green1, ucsf_colors$pink1, ucsf_colors$blue1, ucsf_colors$gray3)

gene_sets_for_sina <- gene_set_data %>% 
  filter(gene_set %in% gene_sets_order) %>% 
  rename('strain' = query)

data <-
  corr_for_sina %>% 
  left_join(gene_sets_for_sina, by = 'strain') %>% 
  filter(greater_fdr < 0.05) %>% 
  mutate(gene_set = ifelse(is.na(gene_set), 'other', gene_set) %>% 
           factor(levels = gene_sets_order)) %>% 
  arrange(mutant_group, gene_set)

sinaplot <- ggplot(data, aes(x = mutant_group, y = pearson, size = greater_fdr)) +
  geom_violin(data = filter(data, gene_set == 'other'),
              fill = ucsf_colors$gray3, color=NA,
              width = 0.8, alpha = 0.3) +
  geom_jitter(data = filter(data, gene_set == 'nuclear pore complex'),
              mapping = aes(x = as.numeric(mutant_group) - 0.3, size = greater_fdr, fill = gene_set),
              width = 0.1, pch = 21) +
  geom_jitter(data = filter(data, gene_set == 'spindle assembly checkpoint'),
              mapping = aes(size = greater_fdr, fill = gene_set), 
              width = 0.1, pch = 21) +
  geom_jitter(data = filter(data, gene_set == 'tRNA modification'),
              mapping = aes(x = as.numeric(mutant_group) + 0.3, size = greater_fdr, fill = gene_set),
              width = 0.1, pch = 21) +
  scale_fill_manual(name='Gene set', values=gene_sets_colors) +
  scale_size_continuous(name = 'P-value', range = c(1.3, 0.4),
                        breaks = c(0.001, 0.05), limits = c(0, 0.06)) +
  ylim(c(0.05, 0.45)) + xlab('Gsp1 mutant group') +
  ylab('Gsp1 mutant - gene\nPearson correlation') +
  guides(size=guide_legend(nrow=4,byrow=TRUE,
                            override.aes = list(fill=NA)),
          color=guide_legend(nrow=3,byrow=TRUE)) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica", size = 6),
    axis.title = element_text(size = 6), axis.text = element_text(size = 6),
    axis.ticks = element_line(size = 0.05), axis.ticks.length = unit(0.05, 'cm'),
    legend.position = 'right',
    legend.direction = 'vertical',
    #legend.position = 'none',
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.01, 'cm'),
    legend.box = 'vertical', 
    legend.box.just = 'left',
    legend.text = element_text(size = 2), 
    legend.title = element_text(size = 2),
    legend.margin = margin(t = 0, unit='cm'),
    plot.margin = margin(t = 0, unit='cm'),
    axis.line = element_line(size = 0.1),
    strip.text.x = element_text(size = 6)
  )
############ SINAPLOT OVER



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

theme_custom = function(gene_set = 'all', show_x_labels = TRUE) {
  theme_classic() %+replace%
    theme(
      axis.text.y = element_text(colour = ifelse(gene_set == 'all', 'white', 'black')),
      axis.text.x = element_text(size = ifelse(show_x_labels, 6, 0.1))
    )
}

### GENERIC FUNCTION TO MAKE HETMAPS PER GENESET
plot_fdr_or_pearson <- function(data, 
                                fill_label = 'FDR of positive Pearson correlation',
                                fill_colors = c(ucsf_colors$gray1, 'white'),
                                fill_limits = c(0, 0.1),
                                fill_breaks = c(0.1, 0),
                                y_axis_label = '',
                                x_axis_label = '',
                                show_x = FALSE,
                                legend_position = 'none'
                                ) {
  data %>% 
    mutate('group' = case_when(mutant %in% mut_groups$I ~ 'I',
                                      mutant %in% mut_groups$II ~ 'II',
                                      mutant %in% mut_groups$III ~ 'III')) %>% 
    ggplot(aes(y = query, x = mutant, fill = value, group = mutant)) + 
    facet_grid(~group, scales = "free_x") +  # scales is important not to repeat all the mutants in each group
    geom_tile(color = ucsf_colors$gray3) + 
    labs(fill = fill_label) +
    scale_fill_gradientn(colors = fill_colors, 
                       limits = fill_limits, breaks = fill_breaks, na.value = 'white') +
    theme_minimal() +
    ylab(y_axis_label) +
    xlab(x_axis_label) +
    theme_custom(gene_set = gene_set_to_plot, show_x_labels = show_x) +
    theme(text = element_text(size = 6),
        axis.title.x = element_text(size = 6, margin = margin(0, 0, 0, 0)),
        axis.title.y = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, margin = margin(0, 0, 0, 0)),
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = 0.05),
        axis.ticks.length = unit(0.05, 'cm'),
        axis.line = element_line(size = 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.position = legend_position,
        legend.box.margin = margin(0,0,0,0)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
         size = guide_legend(title.position = "top", title.hjust = 1))
}

############ MAKE HEATMAPS PER GENE SET
distance_method <- 'euclidean'
main_fig_gene_sets <- c("spindle assembly checkpoint", "nuclear pore complex", "tRNA modification")
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
  if (i == 3) {
    plot[[gene_set_to_plot]][['fdr']] <- to_plot %>% 
      select(mutant, query, 'value' = greater_fdr) %>% 
      plot_fdr_or_pearson(data = ., fill_label = 'FDR of positive Pearson correlation',
                        fill_colors = c(ucsf_colors$gray1, 'white'),
                        fill_limits = c(0, 0.1),
                        fill_breaks = c(0.1, 0),
                        y_axis_label = '',
                        x_axis_label = '',
                        show_x = T)
    plot[[gene_set_to_plot]][['corr']] <- to_plot %>% 
      select(mutant, query, 'value' = pearson) %>% 
      plot_fdr_or_pearson(data = ., fill_label = 'Pearson correlation',
                          fill_colors = c(ucsf_colors$pink1, 'white', ucsf_colors$green1),
                          fill_limits = c(-0.4, 0.4),
                          fill_breaks = c(-0.4, 0, 0.4),
                          y_axis_label = '',
                          x_axis_label = '',
                          show_x = F)
  } else {
    plot[[gene_set_to_plot]][['fdr']] <- to_plot %>% 
      select(mutant, query, 'value' = greater_fdr) %>% 
      plot_fdr_or_pearson(data = ., fill_label = 'FDR of positive Pearson correlation',
                          fill_colors = c(ucsf_colors$gray1, 'white'),
                          fill_limits = c(0, 0.1),
                          fill_breaks = c(0.1, 0),
                          y_axis_label = '',
                          x_axis_label = '',
                          show_x = FALSE)
    plot[[gene_set_to_plot]][['corr']] <- to_plot %>% 
      select(mutant, query, 'value' = pearson) %>% 
      plot_fdr_or_pearson(data = ., fill_label = 'Pearson correlation',
                          fill_colors = c(ucsf_colors$pink1, 'white', ucsf_colors$green1),
                          fill_limits = c(-0.4, 0.4),
                          fill_breaks = c(-0.4, 0, 0.4),
                          y_axis_label = '',
                          x_axis_label = '',
                          show_x = FALSE)
  } 
  heatmap_height <- to_plot %>% pull(query) %>% unique() %>% length() %>%  `/`(7)
  if (heatmap_height < 1) {heatmap_height <- 1} else if (heatmap_height < 3) {heatmap_height <- 1.85}
  plot[[gene_set_to_plot]][['height']] <- heatmap_height
}

### overwrite the height scaling for main figure
plot[["spindle assembly checkpoint"]][['height']] <- 1.9
plot[["nuclear pore complex"]][['height']] <- 1.85
plot[["tRNA modification"]][['height']] <- 1.35

barplot_fixed_height <- 2.5
sinaplot_height <- 0.75
file_height <- sinaplot_height + 1*barplot_fixed_height  + plot[["spindle assembly checkpoint"]][['height']] + plot[["nuclear pore complex"]][['height']] + plot[["tRNA modification"]][['height']]
combn_plots <- plot_grid(sinaplot, 
                         plot[["spindle assembly checkpoint"]][['fdr']],
                         plot[["nuclear pore complex"]][['fdr']],
                         plot[["tRNA modification"]][['fdr']],
                         line_plot,
                         nrow = 5, align = 'v', axis = c('lr'),
                         rel_heights = c(
                           sinaplot_height,
                           plot[["spindle assembly checkpoint"]][['height']]/barplot_fixed_height,
                           plot[["nuclear pore complex"]][['height']]/barplot_fixed_height,
                           plot[["tRNA modification"]][['height']]/barplot_fixed_height,
                           1)
                         )

pdf(str_c(output_directory, '/Fig_4BC_line_plot.pdf'), width = 3.5, height = file_height)
print(combn_plots)
dev.off()

combn_plots <- plot_grid( 
                         plot[["spindle assembly checkpoint"]][['corr']],
                         plot[["nuclear pore complex"]][['corr']],
                         plot[["tRNA modification"]][['corr']],
                         line_plot,
                         nrow = 4, align = 'v', axis = c('lr'),
                         rel_heights = c(
                           plot[["spindle assembly checkpoint"]][['height']]/barplot_fixed_height,
                           plot[["nuclear pore complex"]][['height']]/barplot_fixed_height,
                           plot[["tRNA modification"]][['height']]/barplot_fixed_height,
                           1)
)
pdf(str_c(output_directory, '/Ext_Data_Fig_9b.pdf'), width = 3.55, height = file_height-sinaplot_height)
print(combn_plots)
dev.off()
# 
EDF_gene_sets <- c("Gsp1 partner", "nuclear transport of protein and mRNA",
                        "CVT pathway mitochondria and vacuole", "actin, tubulin, and polarity",
                        "transcription regulation and mediator complex", "5' mRNA capping",
                        'histones and chromatin')
plot <- list()

for (i in seq_along(EDF_gene_sets)) {
  gene_set_to_plot <- EDF_gene_sets[i]
  to_cluster <- corr %>%
    filter(gene_set == gene_set_to_plot) %>% 
    select(mutant, query, pearson)
  queries_ordered <- pca_order_query_genes_by_pearson(data = to_cluster, dist_method = distance_method)
  to_plot <- corr %>% 
    filter(gene_set == gene_set_to_plot) %>% 
    select(mutant, query, greater_fdr) %>% 
    mutate("mutant" = factor(mutant, mutants_ordered), 
           "query" = factor(query, queries_ordered)
           ) %>% 
    arrange(query, mutant)
  if (i %% 2 == 0) {  ### for even numbers of i
    plot[[gene_set_to_plot]][['fdr']] <- to_plot %>% 
      select(mutant, query, 'value' = greater_fdr) %>% 
      plot_fdr_or_pearson(data = ., fill_label = 'FDR of positive Pearson correlation',
                          fill_colors = c(ucsf_colors$gray1, 'white'),
                          fill_limits = c(0, 0.1),
                          fill_breaks = c(0.1, 0),
                          y_axis_label = '',
                          x_axis_label = '',
                          show_x = T)
  } else {
    plot[[gene_set_to_plot]][['fdr']] <- to_plot %>% 
      select(mutant, query, 'value' = greater_fdr) %>% 
      plot_fdr_or_pearson(data = ., fill_label = 'FDR of positive Pearson correlation',
                          fill_colors = c(ucsf_colors$gray1, 'white'),
                          fill_limits = c(0, 0.1),
                          fill_breaks = c(0.1, 0),
                          y_axis_label = '',
                          x_axis_label = '',
                          show_x = FALSE)
  } 
  heatmap_height <- to_plot %>% pull(query) %>% unique() %>% length() %>%  `/`(8)
  if (heatmap_height < 1) {heatmap_height <- 1.3} #else if (heatmap_height < 3) {heatmap_height <- 2}
  plot[[gene_set_to_plot]][['height']] <- heatmap_height
}



file_height <-  barplot_fixed_height + plot[["Gsp1 partner"]][['height']] + plot[["nuclear transport of protein and mRNA"]][['height']]
combn_plots <- plot_grid(plot[["Gsp1 partner"]][['fdr']],
                         plot[["nuclear transport of protein and mRNA"]][['fdr']],
                         line_plot,
                         nrow = 3, align = 'v', axis = c('lr'),
                         rel_heights = c(
                           (plot[["Gsp1 partner"]][['height']])/barplot_fixed_height,
                           plot[["nuclear transport of protein and mRNA"]][['height']]/barplot_fixed_height,
                           1))
pdf(str_c(output_directory,'/Ext_Fig10a.pdf'), width = 3.3, height = file_height)
print(combn_plots)
dev.off()



# 

file_height <- barplot_fixed_height + plot[["transcription regulation and mediator complex"]][['height']] + plot[["5' mRNA capping"]][['height']]
combn_plots <- plot_grid(plot[["transcription regulation and mediator complex"]][['fdr']],
                         plot[["5' mRNA capping"]][['fdr']],
                         line_plot,
                         nrow = 3, align = 'v', axis = c('lr'),
                         rel_heights = c(
                           (plot[["transcription regulation and mediator complex"]][['height']])/barplot_fixed_height,
                           plot[["5' mRNA capping"]][['height']]/barplot_fixed_height,
                           1))
pdf(str_c(output_directory, '/Ext_Fig10b.pdf'), width = 3.3, height = file_height)
print(combn_plots)
dev.off()


file_height <- barplot_fixed_height + plot[["CVT pathway mitochondria and vacuole"]][['height']] + plot[["actin, tubulin, and polarity"]][['height']]
combn_plots <- plot_grid(plot[["CVT pathway mitochondria and vacuole"]][['fdr']],
                         plot[["actin, tubulin, and polarity"]][['fdr']],
                         line_plot,
                         nrow = 3, align = 'v', axis = c('lr'),
                         rel_heights = c(
                           (plot[["CVT pathway mitochondria and vacuole"]][['height']])/barplot_fixed_height,
                           plot[["actin, tubulin, and polarity"]][['height']]/barplot_fixed_height,
                           1))
pdf(str_c(output_directory, '/Ext_Fig10c.pdf'), width = 3.3, height = file_height)
print(combn_plots)
dev.off()



merged_gene_sets <- c(main_fig_gene_sets, EDF_gene_sets)
by_set_corr <- corr %>% 
  filter(gene_set %in% merged_gene_sets) %>% 
  group_by(mutant, gene_set) %>% 
  summarise('mean_pearson' = mean(pearson, na.rm = T), 
            'median_pearson' = median(pearson, na.rm = T),
            'upper_quartile' = quantile(pearson, na.rm= T)[4]) %>% 
  ungroup()
###### to order the gene sets by mean pearson correlation, use the same function 
### that I used above to order the genes by greater_fdr (adjust the column names :))
gene_sets_ordered_by_mean_pearson <- by_set_corr %>% 
  select(mutant, gene_set, mean_pearson) %>% 
  rename('query' = gene_set, 'greater_fdr' = mean_pearson) %>% 
  pca_order_query_genes(data = ., dist_method = distance_method)

to_plot <- by_set_corr %>% 
  select(mutant, "query" = gene_set, 'value' = upper_quartile) %>% 
  mutate("mutant" = factor(mutant, mutants_ordered), 
         "query" = factor(query, gene_sets_ordered_by_mean_pearson)) %>% 
  arrange(query, mutant)
### then plot using the function for plotting
upper_quartile_plot <- plot_fdr_or_pearson(data = to_plot, 
                    fill_label = 'Q3 Pearson corr',
                    fill_colors = c('white', 'white', ucsf_colors$green1),
                    fill_limits = c(-0.4, 0.4),
                    fill_breaks = c(-0.4, -0.1, 0, 0.4),
                    y_axis_label = '',
                    x_axis_label = '',
                    show_x = TRUE)



to_plot <- by_set_corr %>% 
  select(mutant, "query" = gene_set, 'value' = mean_pearson) %>% 
  mutate("mutant" = factor(mutant, mutants_ordered), 
         "query" = factor(query, gene_sets_ordered_by_mean_pearson)) %>% 
  arrange(query, mutant)
### then plot using the function for plotting
mean_pearson_plot <- plot_fdr_or_pearson(data = to_plot, 
                    fill_label = 'mean Pearson corr',
                    fill_colors = c('white', 'white', ucsf_colors$green1),
                    fill_limits = c(-0.3, 0.3),
                    fill_breaks = c(-0.3, -0.1, 0, 0.3),
                    y_axis_label = '',
                    x_axis_label = '',
                    show_x = TRUE,
                    legend_position = 'right')

pdf(str_c(output_directory, '/Pearson_heatmap_per_gene_set.pdf'), width = 6, height = 3)
print(mean_pearson_plot)
dev.off()

