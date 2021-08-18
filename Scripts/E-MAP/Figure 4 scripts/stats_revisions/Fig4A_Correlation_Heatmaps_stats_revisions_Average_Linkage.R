library(tidyverse)
library(cowplot)
library(factoextra)
library(dendextend)
library(ComplexHeatmap)
library(circlize)

source('ucsf_colors.R')

##### LOAD DATA

# load correlation matrix
load('Data/filtered_v6_correlations.RData')
filtered_correlations <- filter(filtered_correlations, !query_uniq2 %in% c('crm1_damp','yrb2_damp'))

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
row.hc <- clustfn(corr_mat, dist_method = 'euclidean', clust_method = 'average')
col.hc <- clustfn(t(corr_mat), dist_method = 'pearson', clust_method = 'average')


##### ROTATION OF THE TREES

# rotate the yeast query gene tree so the density goes from upper left to lower right
row.hc <- rotate(row.hc, rev(row.hc$labels[row.hc$order]))


figure_order <- c('D79S','T34Q','T34E','K101R','D79A','T34G',
                    'T34A','Q147E','R108I','R108L','G80A','Y157A',
                    'H141E','H141R','R108Y','R108Q','R108G','Y148I',
                    'H141I','R112A','R112S','R78K')

col.dend <- as.dendrogram(col.hc)
col.dend <- rotate(col.dend, figure_order)
col.hc <- as.hclust(col.dend)


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


# plotting heatmap
pval_heatmap <- Heatmap(
  corr_mat, show_heatmap_legend = T, name= 'p-value of\npositive\nPearson\ncorrelation',
  col = colorRamp2(c(0, 0.1), c(ucsf_colors$gray1, 'white')),
  width = unit(4.3, 'cm'), gap = unit(1, "mm"),
  
  # columns (mutants)
  column_title = 'Gsp1 point mutant', column_title_gp = gpar(fontsize = 6),
  cluster_columns = as.dendrogram(col.hc), column_split = 5, column_dend_height = unit(6, "mm"),
  column_names_side = 'top', column_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
  
  # mutant cluster annotations
  top_annotation = HeatmapAnnotation(
    which = 'column', height = unit(2, "mm"),
    annot = anno_block(
      # labels = c('I','II','III'),
      labels = c('I','II','III','IV','V'),
      labels_gp = gpar(col = 'white', fontsize = 6, fontfamily='Helvetica'),
      gp = gpar(fill = c('I' = ucsf_colors$navy2,
                         'II' = ucsf_colors$navy2,
                         'III' = ucsf_colors$navy2,
                         'III' = ucsf_colors$navy2,
                         'III' = ucsf_colors$navy2)
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

)

co <- unlist(column_order(pval_heatmap))
dev.off()

add_error_bars_and_stars <- function(annot, val, err) {
  decorate_annotation(annot, {
    # x positions
    err_w = 0.2  # width of error bars
    x = c(seq(1,13), 14.54, seq(16.08,23.08))
    x0 = c(seq(1,13), 14.54, seq(16.08,23.08)) - err_w
    x1 = c(seq(1,13), 14.54, seq(16.08,23.08)) + err_w
    
    # y positions
    y = val[co]
    y0 = val[co] - err[co]
    y1 = val[co] + err[co]
    
    # draw error bars
    grid.segments(x0=x, x1=x, y0=y0, y1=y1, default.units = 'native', gp = gpar(lex=0.5))
    grid.segments(x0=x0, x1=x1, y0=y0, y1=y0, default.units = 'native', gp = gpar(lex=0.5))
    grid.segments(x0=x0, x1=x1, y0=y1, y1=y1, default.units = 'native', gp = gpar(lex=0.5))
    if (annot == 'in vitro\nGAP relative\nefficiency\n(MUT/WT)') {
      x_star = c(17.08, 18.08, 20.08)
      y_star = rep(0.4, 3)
    } else if (annot == 'in vitro\nGEF relative\nefficiency\n(MUT/WT)') {
      x_star = c(13, 17.08, 18.08, 20.08)
      y_star = rep(0.1, 4)
    } else if (annot == 'Ln ratio of\nGAP/GEF rel.\nefficiences') {
      x_star = c(13, 17.08, 18.08, 20.08)
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
pdf('Revisions2/Fig_and_analysis_for_StatReviewer/SuppFig17_4A_Corr_Pvalue_Heatmap_Average.pdf', width = 3.7, height = 4.5)
draw(pval_heatmap)

add_error_bars_and_stars('in vitro\nGAP relative\nefficiency\n(MUT/WT)',
                         GAP_kcat_Km$rel_to_WT, GAP_kcat_Km$se)

add_error_bars_and_stars('in vitro\nGEF relative\nefficiency\n(MUT/WT)',
                         GEF_kcat_Km$rel_to_WT, GEF_kcat_Km$se)

add_error_bars_and_stars('Ln ratio of\nGAP/GEF rel.\nefficiences',
                         rel_effic$ln_rel_to_WT, rel_effic$ln_se)

dev.off()



