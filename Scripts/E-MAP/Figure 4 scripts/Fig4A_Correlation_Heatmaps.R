library(tidyverse)
library(cowplot)
library(factoextra)
library(dendextend)
library(ComplexHeatmap)
library(circlize)

source('ucsf_colors.R')

##### LOAD DATA

# load correlation matrix
# load('Data/filtered_v6_correlations.RData')

load('Data/filtered_correlations_BNF_2.RData')
filtered_correlations <- filtered_correlations_BNF_2



# load kinetics data for plotting as ratio of relative GAP, GEF efficiencies
kinetics <- read_delim('Data/kinetics_data_relative_to_WT.txt', delim='\t', col_types=cols())


##### CLUSTERING

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

# transform the correlations table into a wide format matrix, rows are mutants, columns are queries
corr_mat <-
  filtered_correlations %>%
  filter(grepl('GSP1', query_uniq1) & (! grepl('GSP1', query_uniq2))) %>%
  mutate('mutant' = substring(query_uniq1, first = 8)) %>%
  select(mutant, 'query' = query_uniq2, greater_fdr) %>% 
  spread(query, greater_fdr) %>%
  column_to_rownames('mutant') %>%
  as.matrix() %>% 
  t()

# perform clustering
row.hc <- clustfn(corr_mat, dist_method = 'euclidean', clust_method = 'ward.D2')
col.hc <- clustfn(t(corr_mat), dist_method = 'pearson', clust_method = 'ward.D2')
dim(corr_mat)

##### ROTATION OF THE TREES

# rotate the yeast query gene tree so the density goes from upper left to lower right
# row.hc <- rotate(row.hc, rev(row.hc$labels[row.hc$order]))

# rotate mutants based on kinetics ordering (ratio of GAP/GEF relative efficiencies)
kinetics_ordering <-
  kinetics %>%
  filter(measure == 'GAP/GEF kcat/Km') %>%
  arrange(rel_to_WT) %>%
  pull(mutant)

mutant_ordering <- colnames(corr_mat)
mutant_ordering <- mutant_ordering[order(match(mutant_ordering, kinetics_ordering))]
col.hc <- rotate(col.hc, mutant_ordering)



# Below, we make additional valid rotations for the purposes of visualization
col.dend <- as.dendrogram(col.hc)

# Fix the rotation of the GAP mutants: the kinetic_ordering rotation doesn't rotate the whole
# first subtree correctly, where the D79A/T34G branch should come after the other GAP mutants
# based on kinetics (perhaps this is because T34E is close to but right after D79A and T34G
# in kinetics, and K101R comes after them in the strict kinetics ordering. We rotate to a valid
# tree with the D79A/T34G branch coming after, as it makes the most sense to have D79S and T34Q
# before D79A and T34G)
subtree <- col.dend[[1]][[1]]
lbls <- labels(subtree)
col.dend[[1]][[1]] <- rotate(subtree, c(lbls[! lbls %in% c('D79A','T34G')], 'D79A', 'T34G'))

# Rotate G80A to be after T34A
subtree <- col.dend[[1]][[2]][[1]]
lbls <- labels(subtree)
col.dend[[1]][[2]][[1]] <- rotate(subtree, c(lbls[lbls != 'G80A'], 'G80A'))

# Rotate R78K to last (few correlations), and rotate the R112S/R112A branch to be next to R78K,
# since they have matching kinetics
lbls <- labels(col.dend)
col.dend <- rotate(col.dend, c(lbls[! lbls %in% c('R112A', 'R112S', 'R78K')], 'R112A', 'R112S', 'R78K'))

# finished with rotations
col.hc <- as.hclust(col.dend)

# save the mutant ordering
mutant_ordering_clustered <- col.hc$labels[col.hc$order]
# write(mutant_ordering_clustered, 'Figure4_Multispecificity/Plots/4B_order_of_mutants.txt')
# write(mutant_ordering_clustered, 'Data/4B_order_of_mutants.txt')


##### PREPARE DATA STRUCTURES FOR ANNOTATIONS

# KINETICS
# prepare dataframes of kinetic values with errors for barplots
# The ordering should be the same as the matrix columns, because the heatmap
# reorders the bars based on the dendrogram we provide
get_ordered_kinetics_df <- function(m) {
  col.hc$labels %>% 
    enframe(value = 'mutant', name = NULL) %>% 
    left_join(kinetics, by = 'mutant', ) %>%
    complete(measure, nesting(mutant), fill = list(rel_to_WT = 0, se = 0, ln_rel_to_WT = 0, ln_se = 0)) %>% 
    filter(measure == m) %>%
    mutate(mutant = factor(mutant, levels = col.hc$labels)) %>%
    arrange(mutant)
}

GAP_kcat_Km <- get_ordered_kinetics_df('GAP_kcat_Km')
GEF_kcat_Km <- get_ordered_kinetics_df('GEF_kcat_Km')
rel_effic <- get_ordered_kinetics_df('GAP/GEF kcat/Km')

# set color parameters for the kinetic barplots
rc_col_fn <- colorRamp2(c(-2, 0, 4.5), c(ucsf_colors$orange1, 'white', ucsf_colors$cyan1))
rc_colors <- rc_col_fn(rel_effic$ln_rel_to_WT)


# GENE SETS
# get gene set information to highlight the strains on the heatmap
gene_sets_to_show <-
  read_delim('Supplementary_Data_Tables/Excel_files/gene_sets_final.txt',
             delim = '\t', col_types = cols()) %>%
  select(query, gene_set) %>%
  rename('strain' = query) %>%
  filter(gene_set %in% c('nuclear pore complex',
                         'spindle assembly checkpoint',
                         'tRNA modification'))

gene_set_df <-
  rownames(corr_mat) %>% 
  enframe(name = NULL, value = 'strain') %>% 
  left_join(gene_sets_to_show, by = 'strain') %>% 
  rownames_to_column('position') %>% 
  mutate(position = as.integer(position)) %>% 
  filter(!is.na(gene_set)) %>% 
  mutate('color' = case_when(gene_set == 'nuclear pore complex' ~ ucsf_colors$green1,
                             gene_set == 'spindle assembly checkpoint' ~ ucsf_colors$pink1,
                             gene_set == 'tRNA modification' ~ ucsf_colors$blue1))



# plotting heatmap
pval_heatmap <- Heatmap(
  corr_mat, show_heatmap_legend = F,
  col = colorRamp2(c(0, 0.1), c(ucsf_colors$gray1, 'white')),
  width = unit(4.3, 'cm'), gap = unit(1, "mm"),
  
  # columns (mutants)
  column_title = 'Gsp1 point mutant', column_title_gp = gpar(fontsize = 6),
  cluster_columns = as.dendrogram(col.hc), column_split = 3, column_dend_height = unit(6, "mm"),
  column_names_side = 'top', column_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
  
  # mutant cluster annotations
  top_annotation = HeatmapAnnotation(
    which = 'column', height = unit(2, "mm"),
    annot = anno_block(
      labels = c('I','II','III'),
      labels_gp = gpar(col = 'white', fontsize = 6, fontfamily='Helvetica'),
      gp = gpar(fill = c('I' = ucsf_colors$navy1,
                         'II' = ucsf_colors$navy2,
                         'III' = ucsf_colors$navy3)
      ),
    )
  ),
  
  # mutant barplot annotations
  bottom_annotation = HeatmapAnnotation(
    which = 'col',
    height = unit(4.1, "cm"), 
    gap = unit(c(4, 4), "mm"),
    annotation_height = unit(c(1, 1, 1.3), "cm"),
    
    `in vitro\nGAP relative\nefficiency\n(MUT/WT)` = anno_barplot(
      x = GAP_kcat_Km$rel_to_WT,
      ylim = c(0, 4),
      gp = gpar(fill = ucsf_colors$gray1, lex=0.5),
      axis_param = list(gp = gpar(fontsize = 6, fontfamily='Helvetica'),
                        at = c(0, 1, 2, 3, 4)
                        )
      ),
    
    `in vitro\nGEF relative\nefficiency\n(MUT/WT)` = anno_barplot(
      x = GEF_kcat_Km$rel_to_WT,
      ylim = c(0, 1.1),
      gp = gpar(fill = ucsf_colors$gray1, lex=0.5),
      axis_param = list(gp = gpar(fontsize = 6, fontfamily='Helvetica'),
                        at = c(0, 0.25, 0.5, 0.75, 1)
                        )
      ),
    `Ln ratio of\nGAP/GEF rel.\nefficiences` = anno_barplot(
      x = rel_effic$ln_rel_to_WT,
      ylim = c(-2.2,6.5),
      gp = gpar(fill = rc_colors, lex=0.5),
      axis_param = list(gp = gpar(fontsize = 6, fontfamily='Helvetica'))),
    
    annotation_name_gp = gpar(fontsize = 6, fontfamily='Helvetica', lineheight = 0.9),
    annotation_name_side = 'left',
    annotation_name_rot = 0
  ),
  
  # rows (genes)
  row_title = 'Saccharomyces cerevisiae genes', row_title_gp = gpar(fontsize = 6),
  cluster_rows = as.dendrogram(row.hc), row_split = 7, row_dend_width = unit(6, "mm"),
  row_order = row.hc$order,
  show_row_names = F,
  # row_labels = label_only_gene_set_strains,
  # row_names_gp = gpar(fontsize = 2, fontfamily = 'Helvetica',
  #                     col = gene_set_label_colors),
  
  # yeast gene cluster annotation
  left_annotation = HeatmapAnnotation(
    which = 'row', width = unit(2, "mm"),
    annot = anno_block(
      labels = c('1','2','3','4','5','6','7'),
      labels_gp = gpar(col = 'white', fontsize = 6, fontfamily='Helvetica'),
      labels_rot = 0,
      gp = gpar(fill = c('1' = 'black','2' = 'black', '3' = 'black', '4' = 'black',
                         '5' = 'black','6' = 'black', '7' = 'black')
      )
    )
  ),
  
  
  # dots to annotate the genes from the gene sets
  right_annotation = HeatmapAnnotation(
    which = 'row', width = unit(1, 'mm'),
    npc = anno_mark(at = gene_set_df$position,
                    labels = gene_set_df$strain,
                    labels_gp = gpar(col = gene_set_df$color, fontsize=0),
                    link_gp = gpar(col = gene_set_df$color, lex = 0.5),
                    link_width = unit(2, 'mm')
                    ),
    show_annotation_name = FALSE
  )
)

co <- unlist(column_order(pval_heatmap))
dev.off()

add_error_bars_and_stars <- function(annot, val, err) {
  decorate_annotation(annot, {
    # x positions
    err_w = 0.2  # width of error bars
    x = c(seq(1,6), seq(7.54, 14.54), seq(16.08,23.08))
    x0 = c(seq(1,6), seq(7.54, 14.54), seq(16.08,23.08)) - err_w
    x1 = c(seq(1,6), seq(7.54, 14.54), seq(16.08,23.08)) + err_w
    
    # y positions
    y = val[co]
    y0 = val[co] - err[co]
    y1 = val[co] + err[co]
    
    # draw error bars
    grid.segments(x0=x, x1=x, y0=y0, y1=y1, default.units = 'native', gp = gpar(lex=0.5))
    grid.segments(x0=x0, x1=x1, y0=y0, y1=y0, default.units = 'native', gp = gpar(lex=0.5))
    grid.segments(x0=x0, x1=x1, y0=y1, y1=y1, default.units = 'native', gp = gpar(lex=0.5))
    if (annot == 'in vitro\nGAP relative\nefficiency\n(MUT/WT)') {
      x_star = c(20.08, 21.08, 22.08)
      y_star = rep(0.4, 3)
    } else if (annot == 'in vitro\nGEF relative\nefficiency\n(MUT/WT)') {
      x_star = c(13.54, 20.08, 21.08, 22.08)
      y_star = rep(0.1, 4)
    } else if (annot == 'Ln ratio of\nGAP/GEF rel.\nefficiences') {
      x_star = c(13.54, 20.08, 21.08, 22.08)
      y_star = rep(0.6, 4)
    }
    
    # draw star
    grid.points(x=x_star, y=y_star, pch=8, size = unit(0.75, 'mm'), gp = gpar(lex=0.5))
    
    # draw line at y = 1 for the GAP plot and the GEF plot

    
    if (annot != 'Ln ratio of\nGAP/GEF rel.\nefficiences') {
      y_dots = 1
      grid.lines(x =  c(0.5, 6.5), c(y_dots, y_dots), gp = gpar(lty = 'dotted'), default.units = "native")
      grid.lines(x =  c(7.1, 14.1), c(y_dots, y_dots), gp = gpar(lty = 'dotted'), default.units = "native")
      grid.lines(x =  c(14.8, 23.4), c(y_dots, y_dots), gp = gpar(lty = 'dotted'), default.units = "native")
    }
    
    
  }
  ) 
}

# save heatmap
pdf('Revisions2/4A_Corr_Pvalue_Heatmap_0629.pdf', width = 4.5, height = 6.7)
draw(pval_heatmap)

add_error_bars_and_stars('in vitro\nGAP relative\nefficiency\n(MUT/WT)',
                         GAP_kcat_Km$rel_to_WT, GAP_kcat_Km$se)

add_error_bars_and_stars('in vitro\nGEF relative\nefficiency\n(MUT/WT)',
                         GEF_kcat_Km$rel_to_WT, GEF_kcat_Km$se)

add_error_bars_and_stars('Ln ratio of\nGAP/GEF rel.\nefficiences',
                         rel_effic$ln_rel_to_WT, rel_effic$ln_se)

dev.off()

# # make legends
# pdf('Figure4_Multispecificity/Plots/4A_Corr_Pvalue_Legend.pdf', width = 1, height = 1)
# Legend(title = 'P-value of pos. corr',
#        col_fun = colorRamp2(c(0, 0.1), c(ucsf_colors$gray1, 'white')),
#        at = c(0, 0.05, 0.1),
#        direction = 'horizontal',
#        grid_height = unit(1, 'mm'),
#        grid_width = unit(0.5, "mm"),
#        title_gp = gpar(fontsize = 6),
#        title_position = 'topcenter',
#        labels_gp = gpar(fontsize = 6)) %>% 
#   draw()
# dev.off()
# 
# pdf('Figure4_Multispecificity/Plots/4B_Kinetics_Legend.pdf', width = 1, height = 1)
# Legend(title = 'log(GAP eff - GEF eff)',
#        col_fun = rc_col_fn,
#        at = c(-2, 0, 3.5),
#        direction = 'horizontal',
#        grid_height = unit(1, 'mm'),
#        grid_width = unit(0.5, "mm"),
#        title_gp = gpar(fontsize = 6),
#        title_position = 'topcenter',
#        labels_gp = gpar(fontsize = 6)) %>% 
#   draw()
# dev.off()


# UNCOMMENT TO SAVE THE GENE SETS, NUMBERED AS THEY ARE IN THE HEATMAP
# clustered_row_dend <- row_dend(pval_heatmap)
# labels(clustered_row_dend)
# 
# clustered_row_dend %>%
#   lapply(labels) %>%
#   enframe() %>%
#   unnest(value) %>%
#   rename('cluster' = name, 'strain' = value) %>%
#   write_csv('Supplementary_Data_Tables/Excel_files/corr_clustering_heatmap_cluster_gene_sets.csv')
# 
# labels(clustered_row_dend[[5]])



# 
# # also plot pcc_heatmap
# pcc_corr_mat <-
#   filtered_correlations %>%
#   filter(grepl('GSP1', query_uniq1) & (! grepl('GSP1', query_uniq2))) %>%
#   mutate('mutant' = substring(query_uniq1, first = 8)) %>%
#   select(mutant, 'query' = query_uniq2, pearson) %>%
#   spread(query, pearson) %>%
#   column_to_rownames('mutant') %>%
#   as.matrix()
# 
# 
# # plotting heatmap
# pcc_heatmap <-
#   Heatmap(
#     pcc_corr_mat, name = 'Pearson corr\ncoefficient', show_heatmap_legend = F,
#     col = colorRamp2(c(-0.4, 0, 0.4), c(ucsf_colors$pink1, 'white', ucsf_colors$green1)),
# 
#     # rows
#     row_title = 'Gsp1 point mutant', row_title_gp = gpar(fontsize = 6),
#     cluster_rows = as.dendrogram(row.hc), row_split = 3, row_dend_width = unit(6, "mm"),
#     row_names_side = 'left', row_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
# 
#     # mutant cluster annotations
#     left_annotation = HeatmapAnnotation(
#       which = 'row',
#       width = unit(2, "mm"),
#       annot = anno_block(
#         gp = gpar(fill = c('1' = ucsf_colors$navy1,
#                            '2' = ucsf_colors$navy2,
#                            '3' = ucsf_colors$navy3)),
#         labels = c('1','2','3'),
#         labels_gp = gpar(col = 'white', fontsize = 6, fontfamily='Helvetica')
#       )
#     ),
# 
#     # mutant barplot annotations
#     right_annotation = HeatmapAnnotation(
#       which = 'row',
#       width = unit(4.5, "cm"),
#       `Log ratio of\nGAP/GEF rel.\nefficiences` = anno_barplot(
#         logGAPGEF, gp = gpar(limits = c(-2, 3.5), fill = GAPGEF_colors),
#         axis_param = list(gp = gpar(fontsize = 6, fontfamily='Helvetica'),
#                           labels_rot = 0)),
#       `in vitro\nGAP relative\nefficiency\nMUT/WT` = anno_barplot(
#         GAP_for_barchart, gp = gpar(limits = c(0, 5), fill = ucsf_colors$gray1),
#         axis_param = list(gp = gpar(fontsize = 6, fontfamily='Helvetica'),
#                           labels_rot = 0)),
#       `in vitro\nGEF relative\nefficiency\nMUT/WT` = anno_barplot(
#         GEF_for_barchart, gp = gpar(limits = c(0, 5), fill = ucsf_colors$gray1),
#         axis_param = list(gp = gpar(fontsize = 6, fontfamily='Helvetica'),
#                           labels_rot = 0)),
#       annotation_name_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
#       annotation_name_side = 'top',
#       annotation_name_rot = 0
#     ),
# 
#     # columns
#     column_title = 'Saccharomyces cerevisiae genes', column_title_gp = gpar(fontsize = 6),
#     cluster_columns = as.dendrogram(col.hc), column_split = 7, column_dend_height = unit(6, "mm"),
#     # cluster_columns = F, column_order = col_pcoor_ordering,
#     column_order = col.hc$order,
#     show_column_names = F,
#     column_names_gp = gpar(fontsize = 2, fontfamily='Helvetica'),
#     
#     # yeast gene cluster annotation
#     top_annotation = HeatmapAnnotation(
#       which = 'column', height = unit(2, "mm"),
#       annot = anno_block(
#         gp = gpar(fill = c('1' = 'black','2' = 'black', '3' = 'black', '4' = 'black',
#                            '5' = 'black','6' = 'black', '7' = 'black')),
#         labels = c('1','2','3','4','5','6','7'),
#         labels_gp = gpar(col = 'white', fontsize = 6, fontfamily='Helvetica')
#       )
#     )
#   )
# 
# # make legends
# pcc_heatmap_legend <-
#   Legend(title = 'Pearson corr. coefficient',
#          col_fun = colorRamp2(c(-0.4, 0, 0.4), c(ucsf_colors$pink1, 'white', ucsf_colors$green1)),
#          at = c(-0.4, 0, 0.4),
#          direction = 'horizontal',
#          grid_height = unit(1, 'mm'),
#          grid_width = unit(0.5, "mm"),
#          title_gp = gpar(fontsize = 6),
#          title_position = 'topcenter',
#          labels_gp = gpar(fontsize = 6))
# 
# # save
# pdf('Extended_Figures/Ext_Fig8_PCC_Heatmap.pdf', width = 7.2, height = 2.7)
# draw(pcc_heatmap)
# dev.off()
# 
# pdf('Extended_Figures/Ext_Fig8_PCC_Legend.pdf', width = 1, height = 1)
# draw(pcc_heatmap_legend)
# dev.off()
# 
# pdf('Extended_Figures/Ext_Fig8_Kinetics_Legend.pdf', width = 1, height = 1)
# draw(barchart_legend)
# dev.off()
