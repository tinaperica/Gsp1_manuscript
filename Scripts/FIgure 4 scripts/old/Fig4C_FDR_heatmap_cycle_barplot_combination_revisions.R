library(tidyverse)
library(factoextra)
library(cowplot)
library(dendextend)
library(scales)
#library(ggnewscale)
source('ucsf_colors.R')
### output directory
output_directory <- 'Revisions/Main Figures/Figure4/heatmap_and_barplot_output'
### input files
load('Data/filtered_v5_correlations.RData')
strong_mutants <- filtered_correlations %>% filter(grepl('GSP1', query_uniq1)) %>% pull(query_uniq1) %>% unique()
mutants_ordered <- read_tsv('Data/order_of_mutants.txt', col_names = F)$X1
gene_set_data <- read_tsv('Data/gene_sets.txt') %>% 
  select(query, gene_name, gene_set)

all_queries <- read_tsv('Data/all_queries.txt') %>% 
  select(query, 'gene_name' = name) %>% 
  mutate('gene_set' = 'all')
gene_set_data <- bind_rows(gene_set_data, all_queries)
gene_sets <- gene_set_data %>% pull(gene_set) %>% unique()

# Note for Tina: This is no longer capped
kinetics <- read_tsv("Data/kinetics_data_relative_to_WT.txt") %>% 
  filter(! measure %in% c('int', 'NMR'))

# only plot mutants that have both GAP and GEF kinetic measurements
mutants_with_complete_kinetics <-
  kinetics %>%
  select(mutant, measure, rel_to_WT) %>% 
  pivot_wider(values_from = rel_to_WT, names_from = 'measure') %>% 
  drop_na(GAP_kcat_Km, GEF_kcat_Km) %>% 
  pull(mutant) %>%
  unique()

# Note for Tina: This is the ordering for plots
mutants_for_figures <- intersect(mutants_ordered, mutants_with_complete_kinetics)

corr <- filtered_correlations %>% 
  filter(query_uniq1 %in% strong_mutants) %>% 
  filter(query_uniq2 %in% unique(gene_set_data$query)) %>% 
  mutate('query_uniq1' = substring(query_uniq1, first = 8)) %>% 
  filter(query_uniq1 %in% mutants_for_figures) %>% 
  inner_join(., gene_set_data, by = c('query_uniq2' = 'query')) %>% 
  select('mutant' = query_uniq1, 'query' = query_uniq2, gene_name, gene_set, greater_fdr, pearson)


kinetics_plotting <- function(data, 
                              ylabel = "", 
                              colorbar_label = "", 
                              scale_fill_colors = c(ucsf_colors$orange1,"white",ucsf_colors$cyan1), 
                              scale_fill_values = c(-2,0,5), 
                              scale_fill_limits = c(-2,5)
                              ) {
  data %>% 
  ggplot(aes(mutant, value, fill = value)) + 
  geom_bar(stat = "identity", color = 'black', width = 0.7, size = 0.2) + 
  scale_fill_gradientn(colours = scale_fill_colors, 
                       values = rescale(scale_fill_values),
                       limits = scale_fill_limits) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 0.5, size = 0.2) + 
  #coord_flip() +
  ylab(ylabel) +
  labs(fill = colorbar_label) +
  theme_light() +
  xlab(element_blank()) +
  theme(text = element_text(size = 6), 
        axis.title.x = element_text(size = 6, margin = margin(0, 0, 0, 0)), 
        axis.text.y = element_text(size = 6),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = 0.05),
        axis.ticks.length = unit(0.05, 'cm'),
        axis.line = element_line(size = 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.position = 'right',
        legend.box.margin = margin(0,0,0,0),
        legend.text = element_text(size = 3)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
         size = guide_legend(title.position = "top", title.hjust = 1))
}
###### MAKE ALL THE KINETICS PLOTS
####### ln(GAP/GEF) kcat/Km
kcat_Km_ratio_to_plot <- kinetics %>% 
  filter(mutant %in% mutants_ordered & measure == 'GAP/GEF kcat/Km') %>% 
  mutate("mutant" = factor(mutant, mutants_ordered)) %>% 
  arrange(mutant) %>% 
  select(mutant, 'value' = ln_rel_to_WT, 'se' = ln_se)
ratio_kcat_Km_plot <- kinetics_plotting(data = kcat_Km_ratio_to_plot, 
                                        ylabel = "ln(relative GAP/GEF efficiency)", 
                                        colorbar_label = "ln(GAP kcat/Km(MUT/WT) / GEF kcat/Km(MUT/WT) )", 
                                        scale_fill_colors = c(ucsf_colors$orange1,"white",ucsf_colors$cyan1), 
                                        scale_fill_values = c(-2,0,5), 
                                        scale_fill_limits = c(-2,5))
### ln(GEF kcat/Km MUT/WT)
GEF_kcat_Km_to_plot <- kinetics %>% 
  filter(mutant %in% mutants_ordered & measure == 'GEF_kcat_Km') %>% 
  mutate("mutant" = factor(mutant, mutants_ordered)) %>% 
  arrange(mutant) %>% 
  select(mutant, 'value' = ln_rel_to_WT, 'se' = ln_se)
GEF_kcat_Km_plot <- kinetics_plotting(data = GEF_kcat_Km_to_plot, 
                                        ylabel = "ln(relative GEF efficiency)", 
                                        colorbar_label = "ln(GEF kcat/Km(MUT/WT) )", 
                                        scale_fill_colors = c("white",ucsf_colors$cyan1), 
                                        scale_fill_values = c(-6,0), 
                                        scale_fill_limits = c(-6,0))

### ln(GEF kcat MUT/WT)
GEF_kcat_to_plot <- kinetics %>% 
  filter(mutant %in% mutants_ordered & measure == 'GEF_kcat') %>% 
  mutate("mutant" = factor(mutant, mutants_ordered)) %>% 
  arrange(mutant) %>% 
  select(mutant, 'value' = ln_rel_to_WT, 'se' = ln_se)
GEF_kcat_plot <- kinetics_plotting(data = GEF_kcat_to_plot, 
                                      ylabel = "ln(relative GEF kcat)", 
                                      colorbar_label = "ln(GEF kcat (MUT/WT) )", 
                                      scale_fill_colors = c("white",ucsf_colors$cyan1), 
                                      scale_fill_values = c(-1.7,1.1), 
                                      scale_fill_limits = c(-1.7,1.1))

### ln(GEF Km MUT/WT)
GEF_Km_to_plot <- kinetics %>% 
  filter(mutant %in% mutants_ordered & measure == 'GEF_Km') %>% 
  mutate("mutant" = factor(mutant, mutants_ordered)) %>% 
  arrange(mutant) %>% 
  select(mutant, 'value' = ln_rel_to_WT, 'se' = ln_se)
GEF_Km_plot <- kinetics_plotting(data = GEF_Km_to_plot, 
                                   ylabel = "ln(relative GEF Km)", 
                                   colorbar_label = "ln(GEF Km (MUT/WT) )", 
                                   scale_fill_colors = c(ucsf_colors$cyan1, "white"), 
                                   scale_fill_values = c(-0.7,6), 
                                   scale_fill_limits = c(-0.7,6))
### ln(GAP kcat/Km MUT/WT)
GAP_kcat_Km_to_plot <- kinetics %>% 
  filter(mutant %in% mutants_ordered & measure == 'GAP_kcat_Km') %>% 
  mutate("mutant" = factor(mutant, mutants_ordered)) %>% 
  arrange(mutant) %>% 
  select(mutant, 'value' = ln_rel_to_WT, 'se' = ln_se)
GAP_kcat_Km_plot <- kinetics_plotting(data = GAP_kcat_Km_to_plot, 
                                      ylabel = "ln(relative GAP efficiency)", 
                                      colorbar_label = "ln(GAP kcat/Km(MUT/WT) )", 
                                      scale_fill_colors = c("white",ucsf_colors$orange1), 
                                      scale_fill_values = c(-3,1.5), 
                                      scale_fill_limits = c(-3,1.5))

### ln(GAP kcat MUT/WT)
GAP_kcat_to_plot <- kinetics %>% 
  filter(mutant %in% mutants_ordered & measure == 'GAP_kcat') %>% 
  mutate("mutant" = factor(mutant, mutants_ordered)) %>% 
  arrange(mutant) %>% 
  select(mutant, 'value' = ln_rel_to_WT, 'se' = ln_se)
GAP_kcat_plot <- kinetics_plotting(data = GAP_kcat_to_plot, 
                                      ylabel = "ln(relative GAP kcat)", 
                                      colorbar_label = "ln(GAP kcat(MUT/WT) )", 
                                      scale_fill_colors = c("white",ucsf_colors$orange1), 
                                      scale_fill_values = c(-1,.5), 
                                      scale_fill_limits = c(-1,.5))

### ln(GAP Km MUT/WT)
GAP_Km_to_plot <- kinetics %>% 
  filter(mutant %in% mutants_ordered & measure == 'GAP_Km') %>% 
  mutate("mutant" = factor(mutant, mutants_ordered)) %>% 
  arrange(mutant) %>% 
  select(mutant, 'value' = ln_rel_to_WT, 'se' = ln_se)
GAP_Km_plot <- kinetics_plotting(data = GAP_Km_to_plot, 
                                   ylabel = "ln(relative GAP Km)", 
                                   colorbar_label = "ln(GAP Km(MUT/WT) )", 
                                   scale_fill_colors = c(ucsf_colors$orange1, "white"), 
                                   scale_fill_values = c(-2,2), 
                                   scale_fill_limits = c(-2,2))


#### HEATMAP FOR MAIN FIGURE
hm_data <- kinetics %>% 
  filter(mutant %in% mutants_for_figures,
         measure %in% c('GAP_kcat_Km', 'GEF_kcat_Km')) %>% 
  select(mutant, measure, ln_rel_to_WT) %>% 
  mutate(mutant = factor(mutant, levels = mutants_for_figures),
         label = case_when(ln_rel_to_WT < 0 ~ '-',
                           ln_rel_to_WT > 0 ~ '+'))

ggplot() +
  geom_tile(data = filter(hm_data, measure == 'GAP_kcat_Km'),
            aes(x = mutant, y = measure, fill = ln_rel_to_WT)) +
  geom_text(data = filter(hm_data, measure == 'GAP_kcat_Km'),
            aes(x = mutant, y = measure, label = label)) +
  scale_fill_gradient2('GAP ln(kcat/Km), MUT/WT',
                       low = ucsf_colors$orange1, mid = 'white', high = ucsf_colors$purple1,
                       guide = guide_colorbar(order = 1)) +
  # scale_fill_gradient('GAP ln(kcat/Km), MUT/WT', low = 'white', high = ucsf_colors$orange1) +
  new_scale_fill() +
  geom_tile(data = filter(hm_data, measure == 'GEF_kcat_Km'),
            aes(x = mutant, y = measure, fill = ln_rel_to_WT)) +
  scale_fill_gradient('GEF ln(kcat/Km), MUT/WT', limits = c(-4, 0),
                      oob=squish,  # very handy, oob = squish sets out-of-bounds vals to limits
                      high = 'white', low = ucsf_colors$cyan1,
                      guide = guide_colorbar(order = 2)) +
  scale_y_discrete(limits=c('GEF_kcat_Km', 'GAP_kcat_Km')) +
  theme(axis.text.x = element_text(angle = 90))

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

plot_fdr_or_pearson <- function(data, 
                                fill_label = 'FDR of positive Pearson correlation',
                                fill_colors = c(ucsf_colors$gray1, 'white'),
                                fill_limits = c(0, 0.1),
                                fill_breaks = c(0.1, 0),
                                y_axis_label = '',
                                x_axis_label = '',
                                show_x = FALSE
                                ) {
  data %>% 
    ggplot(aes(y = query, x = mutant, fill = value)) + 
    geom_tile(color = ucsf_colors$gray3) + 
    labs(fill = fill_label) +
    scale_fill_gradientn(colors = fill_colors, 
                       limits = fill_limits, breaks = fill_breaks, na.value = 'white') +
    theme_minimal() +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "right") +
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
        legend.position = 'none',
        legend.box.margin = margin(0,0,0,0)) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 1),
         size = guide_legend(title.position = "top", title.hjust = 1))
}


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
  if (i == 1) {
    plot[[gene_set_to_plot]][['fdr']] <- to_plot %>% 
      select(mutant, query, 'value' = greater_fdr) %>% 
      plot_fdr_or_pearson(data = ., fill_label = 'FDR of positive Pearson correlation',
                        fill_colors = c(ucsf_colors$gray1, 'white'),
                        fill_limits = c(0, 0.1),
                        fill_breaks = c(0.1, 0),
                        y_axis_label = '',
                        x_axis_label = '',
                        show_x = TRUE)
    plot[[gene_set_to_plot]][['corr']] <- to_plot %>% 
      select(mutant, query, 'value' = pearson) %>% 
      plot_fdr_or_pearson(data = ., fill_label = 'Pearson correlation',
                          fill_colors = c(ucsf_colors$pink1, 'white', ucsf_colors$green1),
                          fill_limits = c(-0.4, 0.4),
                          fill_breaks = c(-0.4, 0, 0.4),
                          y_axis_label = '',
                          x_axis_label = '',
                          show_x = TRUE)
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
  # heatmap_height <- to_plot %>% pull(query) %>% unique() %>% length() %>%  `/`(6)
  # if (heatmap_height < 1.5) {heatmap_height <- 1.5} #else if (heatmap_height < 3) {heatmap_height <- 1.7}
  # plot[[gene_set_to_plot]][['height']] <- heatmap_height
}

# barplot_fixed_height <- 1.5
# file_height <- 1*barplot_fixed_height + plot[["nuclear pore complex"]][['height']] + plot[["spindle assembly checkpoint"]][['height']] + plot[["tRNA modification"]][['height']]
# combn_plots <- plot_grid(plot[["nuclear pore complex"]][['fdr']], 
#                          plot[["spindle assembly checkpoint"]][['fdr']], 
#                          plot[["tRNA modification"]][['fdr']], 
#                          GAP_kcat_Km_plot, GEF_kcat_Km_plot, ratio_kcat_Km_plot,
#                          nrow = 6, align = 'v', axis = c('lr'), 
#                          rel_heights = c(
#                            (1.2*plot[["nuclear pore complex"]][['height']])/barplot_fixed_height, 
#                           plot[["spindle assembly checkpoint"]][['height']]/barplot_fixed_height, 
#                           plot[["tRNA modification"]][['height']]/barplot_fixed_height, 
#                           1, 1, 1))
# file_height <- plot[["spindle assembly checkpoint"]][['height']] +
#                plot[["nuclear pore complex"]][['height']] +
#                plot[["tRNA modification"]][['height']]
# combn_plots <- plot_grid(plot[["spindle assembly checkpoint"]][['fdr']], 
#                          plot[["nuclear pore complex"]][['fdr']], 
#                          plot[["tRNA modification"]][['fdr']], 
#                          nrow = 3, align = 'v', axis = c('lr'), 
                         # rel_heights = c(
                         #   plot[["spindle assembly checkpoint"]][['height']],
                         #   (1.2*plot[["nuclear pore complex"]][['height']]),
                         #   plot[["tRNA modification"]][['height']], 
                         #   1, 1, 1))
                         #   

combn_plots <- plot_grid(plot[["spindle assembly checkpoint"]][['fdr']],
                         plot[["nuclear pore complex"]][['fdr']], 
                         plot[["tRNA modification"]][['fdr']],
                         nrow = 3, align = 'v', axis = c('lr'),
                         rel_heights = c(21,17,9))

pdf(str_c(output_directory, '/Main_Fig4c.pdf'), width = 3, height = file_height*0.8)
print(combn_plots)
dev.off()

file_height <- 1*barplot_fixed_height + plot[["nuclear pore complex"]][['height']] + plot[["spindle assembly checkpoint"]][['height']] + plot[["tRNA modification"]][['height']]
combn_plots <- plot_grid(plot[["nuclear pore complex"]][['corr']], 
                         plot[["spindle assembly checkpoint"]][['corr']], 
                         plot[["tRNA modification"]][['corr']], 
                         GAP_kcat_Km_plot, GEF_kcat_Km_plot, ratio_kcat_Km_plot,
                         nrow = 6, align = 'v', axis = c('lr'), 
                         rel_heights = c(
                           (1.2*plot[["nuclear pore complex"]][['height']])/barplot_fixed_height, 
                           plot[["spindle assembly checkpoint"]][['height']]/barplot_fixed_height, 
                           plot[["tRNA modification"]][['height']]/barplot_fixed_height, 
                           1, 1, 1))
pdf(str_c(output_directory, '/Ext_Data_Fig_9b.pdf'), width = 5, height = file_height)
print(combn_plots)
dev.off()

EDF_gene_sets <- c("Gsp1 partner", "nuclear transport of protein and mRNA", 
                        "CVT pathway mitochondria and vacuole", "actin and polarity", 
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
  if (i %in% c(1, 3, 5, 7)) {
    plot[[gene_set_to_plot]][['fdr']] <- to_plot %>% 
      select(mutant, query, 'value' = greater_fdr) %>% 
      plot_fdr_or_pearson(data = ., fill_label = 'FDR of positive Pearson correlation',
                          fill_colors = c(ucsf_colors$gray1, 'white'),
                          fill_limits = c(0, 0.1),
                          fill_breaks = c(0.1, 0),
                          y_axis_label = '',
                          x_axis_label = '',
                          show_x = TRUE)
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
  heatmap_height <- to_plot %>% pull(query) %>% unique() %>% length() %>%  `/`(6)
  if (heatmap_height < 1.5) {heatmap_height <- 1.5} #else if (heatmap_height < 3) {heatmap_height <- 1.7}
  plot[[gene_set_to_plot]][['height']] <- heatmap_height
}


file_height <- 1*barplot_fixed_height + plot[["Gsp1 partner"]][['height']] + plot[["nuclear transport of protein and mRNA"]][['height']]
combn_plots <- plot_grid(plot[["Gsp1 partner"]][['fdr']], 
                         plot[["nuclear transport of protein and mRNA"]][['fdr']],
                         GAP_kcat_Km_plot, GEF_kcat_Km_plot, ratio_kcat_Km_plot,
                         nrow = 5, align = 'v', axis = c('lr'), 
                         rel_heights = c(
                           (1.2*plot[["Gsp1 partner"]][['height']])/barplot_fixed_height, 
                           plot[["nuclear transport of protein and mRNA"]][['height']]/barplot_fixed_height,
                           1, 1, 1))
pdf(str_c(output_directory, '/Ext_Fig10a.pdf'), width = 5, height = file_height)
print(combn_plots)
dev.off()


file_height <- 1*barplot_fixed_height + plot[["CVT pathway mitochondria and vacuole"]][['height']] + plot[["actin and polarity"]][['height']]
combn_plots <- plot_grid(plot[["CVT pathway mitochondria and vacuole"]][['fdr']], 
                         plot[["actin and polarity"]][['fdr']],
                         GAP_kcat_Km_plot, GEF_kcat_Km_plot, ratio_kcat_Km_plot,
                         nrow = 5, align = 'v', axis = c('lr'), 
                         rel_heights = c(
                           (1.2*plot[["CVT pathway mitochondria and vacuole"]][['height']])/barplot_fixed_height, 
                           plot[["actin and polarity"]][['height']]/barplot_fixed_height,
                           1, 1, 1))
pdf(str_c(output_directory, '/Ext_Fig10b.pdf'), width = 5, height = file_height)
print(combn_plots)
dev.off()


file_height <- 1*barplot_fixed_height + plot[["transcription regulation and mediator complex"]][['height']] + plot[["5' mRNA capping"]][['height']]
combn_plots <- plot_grid(plot[["transcription regulation and mediator complex"]][['fdr']], 
                         plot[["5' mRNA capping"]][['fdr']],
                         GAP_kcat_Km_plot, GEF_kcat_Km_plot, ratio_kcat_Km_plot,
                         nrow = 5, align = 'v', axis = c('lr'), 
                         rel_heights = c(
                           (1.2*plot[["transcription regulation and mediator complex"]][['height']])/barplot_fixed_height, 
                           plot[["5' mRNA capping"]][['height']]/barplot_fixed_height,
                           1, 1, 1))
pdf(str_c(output_directory, '/Ext_Fig10c.pdf'), width = 5, height = file_height)
print(combn_plots)
dev.off()

