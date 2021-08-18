library(tidyverse)

library(factoextra)
library(dendextend)
library(ComplexHeatmap)
library(circlize)

library(gridExtra)
library(grid)

source('ucsf_colors.R')

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


make_heatmap <- function(corr_mat, name, title, i) {

  set.seed(3)
  row.hc <- clustfn(corr_mat, dist_method = 'euclidean', clust_method = 'ward.D2')
  col.hc <- clustfn(t(corr_mat), dist_method = 'pearson', clust_method = 'ward.D2')

  row.hc <- rotate(row.hc, rev(row.hc$labels[row.hc$order]))


  figure_order <- c('D79S','T34Q','T34E','K101R','D79A','T34G',
                    'T34A','Q147E','R108I','R108L','G80A','Y157A',
                    'H141E','H141R','R108Y','R108Q','R108G','Y148I',
                    'H141I','R112A','R112S','R78K')

  col.dend <- as.dendrogram(col.hc)
  col.dend <- rotate(col.dend, figure_order)
  col.hc <- as.hclust(col.dend)



  # plotting heatmap
  hm <- Heatmap(
    corr_mat, show_heatmap_legend = F, 
    col = colorRamp2(c(0, 0.1), c(ucsf_colors$gray1, 'white')),
    width = unit(4.3, 'cm'), gap = unit(1, "mm"),
    
    # columns (mutants)
    column_title = paste(name, '', 'Gsp1 point mutant', sep='\n'),
    
    column_title_gp = gpar(fontsize = 6),
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
                          'III' = ucsf_colors$navy3)))),

    # rows (genes)
    row_title = 'S. cerevisiae alleles', row_title_gp = gpar(fontsize = 6),
    cluster_rows = as.dendrogram(row.hc), row_split = 7, row_dend_width = unit(6, "mm"),
    row_order = row.hc$order, show_row_names = F,
  
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
  )
    return(grid.grabExpr(draw(hm)))
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


# for each step i in the loop, remove i * 10% of data and cluster
heatmaps <- list()
for (i in seq(0,5)) {
  n <- floor((1-0.1*i) * nrow(corr_mat))
  rows_subset <- sample(nrow(corr_mat), n, replace=FALSE)
  corr_mat_tmp <- corr_mat[rows_subset,]
  percent <- 100-(10*i)
  name <- paste0(percent,'% of yeast alleles used, n = ', n)
  heatmaps[[i+1]] <- make_heatmap(corr_mat_tmp, name=name, title=title, i=i)
}


# plot the heatmap, with slice
pdf('Revisions2/Fig_and_analysis_for_StatReviewer/SuppFig20_4A_Corrs_data_subsampling.pdf', height = 6, width = 7)
grid.draw(arrangeGrob(grobs=heatmaps, nrow=2))
dev.off()
