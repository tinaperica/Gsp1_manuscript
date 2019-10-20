### This code makes heatmaps based on correlation data (FDR or Pearson heatmaps) for Figure 1

library(tidyverse)
library(dendextend)
library(factoextra)
library(ComplexHeatmap)
source('ucsf_colors.R')
# load('Data/filtered_v5_correlations.RData')
load('Data/filtered_v6_correlations.RData')


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

# transform the table into a wide format matrix, rows are mutants, columns are queries
corr_mat <-
  filtered_correlations %>%
  filter(grepl('GSP1', query_uniq1) & (! grepl('GSP1', query_uniq2))) %>%
  select('mutant' = query_uniq1, 'query' = query_uniq2, everything()) %>%
  mutate('mutant' = substring(mutant, first = 8)) %>%
  select(mutant, query, greater_fdr) %>%
  spread(query, greater_fdr) %>%
  column_to_rownames('mutant') %>%
  as.matrix()

##### ward.D2 clustering
row.hc <- clustfn(corr_mat, dist_method = 'pearson', clust_method = 'ward.D2')
col.hc <- clustfn(t(corr_mat), dist_method = 'pearson', clust_method = 'ward.D2')

# FOR filtered_corr_v5 ALSO USE THIS TREE FLIPPING
# flipping the first subtree means going into order and subtracting
# the position from nleaves(subtree1)+1. So if there are 100 leaves in
# subtree1, then the gene in position 1 should be in position
# 100+1-1 = 100, the gene in position 10 should be in 100+1-10 = 91,
# and so on.
# 
# subtree1 = as.dendrogram(col.hc)[[1]][[1]]
# ordering <- col.hc$labels[col.hc$order]
# ordering[1:nleaves(subtree1)] <- rev(ordering[1:nleaves(subtree1)])
# col.hc <- rotate(col.hc, ordering)

# reverse the whole tree
col.hc <- rotate(col.hc, rev(col.hc$labels[col.hc$order]))

pdf('E-MAP_Figure1/E-MAP_Figure1/1E_Positive_FDR_wardD2.pdf', width = 4.35, height = 2.5)
# pdf('E-MAP_Figure1/E-MAP_Figure1/1E_Positive_FDR_wardD2_v5.pdf', width = 4.5, height = 2.5)
Heatmap(
  corr_mat, name = 'FDR corr', show_heatmap_legend = F,
  col = colorRamp2(c(0, 0.1), c(ucsf_colors$purple1, 'white')),

  # rows
  row_title = 'Gsp1 point mutant', row_title_gp = gpar(fontsize = 6),
  cluster_rows = as.dendrogram(row.hc), row_split = 3, row_dend_width = unit(6, "mm"),
  row_names_side = 'left', row_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
  left_annotation = HeatmapAnnotation(
    which = 'row', width = unit(2, "mm"),
    annot = anno_block(
      gp = gpar(fill = c('1' = ucsf_colors$navy1,
                         '2' = ucsf_colors$navy2,
                         '3' = ucsf_colors$navy3)),
      labels = c('1','2','3'),
      labels_gp = gpar(col = 'white', fontsize = 4, fontfamily='Helvetica')
    )
  ),

  # columns
  column_title = 'Yeast gene', column_title_gp = gpar(fontsize = 6),
  cluster_columns = as.dendrogram(col.hc), column_split = 4, column_dend_height = unit(6, "mm"),
  column_order = col.hc$order,
  show_column_names = F,
  column_names_gp = gpar(fontsize = 2, fontfamily='Helvetica'),
  top_annotation = HeatmapAnnotation(
    which = 'column', height = unit(2, "mm"),
    annot = anno_block(
      gp = gpar(fill = c('1' = 'black','2' = ucsf_colors$gray1,
                         '3' = ucsf_colors$gray2,'4' = ucsf_colors$gray3)),
      labels = c('1','2','3','4'),
      labels_gp = gpar(col = 'white', fontsize = 4, fontfamily='Helvetica')
    )
  )
)
dev.off()


pdf('E-MAP_Figure1/E-MAP_Figure1/1E_legend.pdf', width = 1, height = 1)
draw(Legend(title = 'FDR',
            col_fun = colorRamp2(c(0, 0.1), c(ucsf_colors$purple1, 'white')),
            at = c(0, 0.05, 0.1),
            direction = 'horizontal',
            grid_height = unit(1, 'mm'),
            grid_width = unit(0.5, "mm"),
            title_gp = gpar(fontsize = 6),
            title_position = 'topcenter',
            labels_gp = gpar(fontsize = 6)
))
dev.off()


# ##### single clustering
# n_clust = 4
# group = cutree(col.hc, k = n_clust)
# col.hc <- clustfn(t(corr_mat), dist_method = 'pearson', clust_method = 'single')
# 
# pdf('E-MAP_Figure1/E-MAP_Figure1/1E_Positive_FDR_single.pdf', width = 9, height = 4)
# Heatmap(
#   corr_mat, name = 'FDR corr', show_heatmap_legend = F,
#   col = colorRamp2(c(0, 0.1), c(ucsf_colors$purple1, 'white')),
# 
#   # rows
#   row_title = 'Gsp1 mutant', row_title_gp = gpar(fontsize = 6),
#   cluster_rows = as.dendrogram(row.hc), row_split = 3, row_dend_width = unit(6, "mm"),
#   row_names_side = 'left', row_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
#   left_annotation = HeatmapAnnotation(
#     which = 'row',
#     annot = anno_block(
#       gp = gpar(fill = c('1' = ucsf_colors$orange1,
#                          '2' = ucsf_colors$gray1,
#                          '3' = ucsf_colors$cyan1)),
#       labels = c('1','2','3'),
#       labels_gp = gpar(col = 'white', fontsize = 6, fontfamily='Helvetica')
#     )
#   ),
# 
#   # columns
#   column_title = 'Yeast gene', column_title_gp = gpar(fontsize = 6),
#   cluster_columns = as.dendrogram(col.hc), column_split = 4, column_dend_height = unit(8, "mm"),
#   column_order = col.hc$order, column_names_gp = gpar(fontsize = 2, fontfamily='Helvetica'),
#   top_annotation = HeatmapAnnotation(
#     which = 'column',
#     annot = group, simple_anno_size = unit(2, "mm"),
#     show_annotation_name = F, show_legend = F,
#     col = list(annot = c('1' = ucsf_colors$cyan1,'2' = ucsf_colors$orange1,
#                          '3' = ucsf_colors$navy1,'4' = ucsf_colors$pink1))
#   )
# )
# dev.off()
