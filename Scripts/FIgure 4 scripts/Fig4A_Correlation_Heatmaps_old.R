library(tidyverse)
library(cowplot)
library(factoextra)
library(dendextend)
library(ComplexHeatmap)
library(circlize)

source('ucsf_colors.R')

# load correlation matrix
load('Data/filtered_v6_correlations.RData')

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

# prepare kinetics data for plotting as ratio of relative GAP, GEF efficiencies
GAP_kinetics <- read_tsv('Data/RanGAP_assay/GAP_kinetics_MichaelisMenten_parameters.txt', col_types = cols()) %>% 
  select(mutant, 'kcat_Km' = mean_kcat_Km, 'kcat' = mean_kcat, 'Km' = mean_Km) %>% unique()
GEF_kinetics <- read_tsv('Data/RanGEF_assay/GEF_kinetics_MichaelisMenten_parameters.txt', col_types = cols())

WT_GEF <- filter(GEF_kinetics, mutant == 'WT')
WT_GAP <- filter(GAP_kinetics, mutant == 'WT')

GEF_kinetics <-
  GEF_kinetics %>%
  mutate('rel_GEF_kcat_Km' = kcat_Km/WT_GEF$kcat_Km) %>%
  mutate('rel_GEF_kcat' = kcat/WT_GEF$kcat) %>%
  mutate('rel_GEF_Km' = Km/WT_GEF$Km) %>%
  select(mutant, rel_GEF_kcat_Km, rel_GEF_kcat, rel_GEF_Km)

GAP_kinetics <-
  GAP_kinetics %>%
  mutate('rel_GAP_kcat_Km' = kcat_Km/WT_GAP$kcat_Km) %>%
  mutate('rel_GAP_kcat' = kcat/WT_GAP$kcat) %>%
  mutate('rel_GAP_Km' = Km/WT_GAP$Km) %>%
  select(mutant, rel_GAP_kcat_Km, rel_GAP_kcat, rel_GAP_Km)

kinetics <-
  inner_join(GEF_kinetics, GAP_kinetics, by = 'mutant') %>%
  mutate('GAP/GEF' = rel_GAP_kcat_Km/rel_GEF_kcat_Km) %>%
  mutate('GAP/GEF' = ifelse(`GAP/GEF` > 30, 30, `GAP/GEF`)) %>%
  select(mutant, `GAP/GEF`)

# get mutant ordering based on ratio of GAP GEF relative efficiencies
GAPGEF_ordering <-
  kinetics %>%
  arrange(`GAP/GEF`) %>%
  pull(mutant)

# transform the table into a wide format matrix, rows are mutants, columns are queries
corr_mat <-
  filtered_correlations %>%
  filter(grepl('GSP1', query_uniq1) & (! grepl('GSP1', query_uniq2))) %>%
  mutate('mutant' = substring(query_uniq1, first = 8)) %>%
  select(mutant, 'query' = query_uniq2, greater_fdr) %>% 
  spread(query, greater_fdr) %>%
  column_to_rownames('mutant') %>%
  as.matrix()

corr_mat %>% dim()

# perform clustering
row.hc <- clustfn(corr_mat, dist_method = 'pearson', clust_method = 'ward.D2')
col.hc <- clustfn(t(corr_mat), dist_method = 'euclidean', clust_method = 'ward.D2')

mutant_ordering <- rownames(corr_mat)
mutant_ordering <- mutant_ordering[order(match(mutant_ordering,GAPGEF_ordering))]

# rotate mutants based on GAP GEF kinetics ordering
row.hc <- rotate(row.hc, mutant_ordering)

# Rotate R78K to last (few correlations)
row.hc <- rotate(row.hc, c(mutant_ordering[mutant_ordering != 'R78K'], 'R78K'))

# Rotate G80A to be after T34A
row.dend <- as.dendrogram(row.hc)
subtree <- row.dend[[1]][[2]][[1]]
lbls <- labels(subtree)
row.dend[[1]][[2]][[1]] <- rotate(subtree, c(lbls[lbls != 'G80A'], 'G80A'))
row.hc <- as.hclust(row.dend)

ordered_of_mutants_to_save <- row.hc$labels[row.hc$order]
write(ordered_of_mutants_to_save, 'Figure4_Multispecificity/Plots/4B_order_of_mutants.txt')

# col_pcoor_ordering <-
#   cmdscale(get_dist(t(corr_mat), method = 'euclidean'), eig = T, k = 1)$point %>%
#   as_tibble(rownames = 'var') %>%
#   arrange(`V1`) %>%
#   pull(var) %>% 
#   rev()

# rotate the yeast query gene tree so the density goes from upper left to lower right
col.hc <- rotate(col.hc, rev(col.hc$labels[col.hc$order]))
# col.dend <- as.dendrogram(col.hc)
# subtree <- col.dend[[2]]
# col_subtree_reordering <- c(labels(subtree[[2]]), rev(labels(subtree[[1]][[1]])), labels(subtree[[1]][[2]]))
# col.dend[[2]] <- rotate(subtree, col_subtree_reordering)
# col.hc <- as.hclust(col.dend)

# get the logGAPGEF values for plotting the barchart alongside the heatmap
logGAPGEF <-
  data.frame(row.hc$labels) %>% 
  mutate('mutant' = as.character(row.hc.labels)) %>% 
  left_join(kinetics) %>% 
  mutate(`GAP/GEF` = ifelse(is.na(`GAP/GEF`), 1, `GAP/GEF`)) %>% 
  mutate(logGAPGEF = log(`GAP/GEF`)) %>%
  pull(logGAPGEF)

GAPGEF_col_fn <- colorRamp2(c(-2, 0, 3.5), c(ucsf_colors$orange1, 'white', ucsf_colors$cyan1))
GAPGEF_colors <- GAPGEF_col_fn(logGAPGEF)

# get just GAP and just GEF efficiencies for bar charts
GAP_for_barchart <-
  data.frame(row.hc$labels) %>%
  mutate('mutant' = as.character(row.hc.labels)) %>%
  left_join(GAP_kinetics) %>%
  mutate(`GAP relative efficiency` = rel_GAP_kcat_Km) %>%
  pull(`GAP relative efficiency`)

GEF_for_barchart <-
  data.frame(row.hc$labels) %>%
  mutate('mutant' = as.character(row.hc.labels)) %>%
  left_join(GEF_kinetics) %>%
  mutate(`GEF relative efficiency` = rel_GEF_kcat_Km) %>%
  pull(`GEF relative efficiency`)

# get partner strains for annotating columns of heatmap
name2ORF <- read_tsv('Data/spitzemap_name2ORF_index.txt')
partner_index <- data.frame(
  name = c('MSN5','SRP1','LOS1','YRB1','YRB2','KAP95','RNA1',
           'SRM1','MTR10','PSE1','NTF2','CRM1','CSE1','KAP104'),
  ORF = c('YDR335W','YNL189W','YKL205W','YDR002W','YIL063C','YLR347C','YMR235C',
          'YGL097W','YOR160W','YMR308C','YER009W','YGR218W','YGL238W','YBR017C'))
partner_strains <- filter(name2ORF, str_detect(ORF, paste(partner_index$ORF, collapse = '|'))) %>% pull(name)
                   
# plotting heatmap
pval_heatmap <- Heatmap(
  corr_mat, name = 'P-value of\npositive correlation', show_heatmap_legend = F,
  col = colorRamp2(c(0, 0.1), c(ucsf_colors$purple1, 'white')),
  
  # rows
  row_title = 'Gsp1 point mutant', row_title_gp = gpar(fontsize = 6),
  cluster_rows = as.dendrogram(row.hc), row_split = 3, row_dend_width = unit(6, "mm"),
  row_names_side = 'left', row_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
  
  # mutant cluster annotations
  left_annotation = HeatmapAnnotation(
    which = 'row',
    width = unit(2, "mm"),
    annot = anno_block(
      gp = gpar(fill = c('1' = ucsf_colors$navy1,
                         '2' = ucsf_colors$navy2,
                         '3' = ucsf_colors$navy3)),
      labels = c('1','2','3'),
      labels_gp = gpar(col = 'white', fontsize = 6, fontfamily='Helvetica')
    )
  ),
  
  # mutant barplot annotations
  right_annotation = HeatmapAnnotation(
    which = 'row',
    width = unit(4, "cm"),
    `Log ratio of\nGAP/GEF rel.\nefficiences` = anno_barplot(
      x = logGAPGEF, gp = gpar(limits = c(-2, 3.5), fill = GAPGEF_colors),
      axis_param = list(gp = gpar(fontsize = 6, fontfamily='Helvetica'),
                        labels_rot = 0)),
    `in vitro\nGAP relative\nefficiency\nMUT/WT` = anno_barplot(
      x = GAP_for_barchart, gp = gpar(limits = c(0, 5), fill = ucsf_colors$gray1),
      axis_param = list(gp = gpar(fontsize = 6, fontfamily='Helvetica'),
                        labels_rot = 0)),
    `in vitro\nGEF relative\nefficiency\nMUT/WT` = anno_barplot(
      x = GEF_for_barchart, gp = gpar(limits = c(0, 5), fill = ucsf_colors$gray1),
      axis_param = list(gp = gpar(fontsize = 6, fontfamily='Helvetica'),
                        labels_rot = 0)),
    
    annotation_name_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
    annotation_name_side = 'top',
    annotation_name_rot = 0
  ),
  
  # columns
  column_title = 'Saccharomyces cerevisiae genes', column_title_gp = gpar(fontsize = 6),
  cluster_columns = as.dendrogram(col.hc), column_split = 7, column_dend_height = unit(6, "mm"),
  column_order = col.hc$order,
  # cluster_columns = F, column_order = col_pcoor_ordering,
  # show_column_names = T,
  column_names_gp = gpar(fontsize = 2, fontfamily = 'Helvetica'),
  
  # yeast gene cluster annotation
  top_annotation = HeatmapAnnotation(
    which = 'column', height = unit(2, "mm"),
    annot = anno_block(
      gp = gpar(fill = c('1' = 'black','2' = 'black', '3' = 'black', '4' = 'black',
                         '5' = 'black','6' = 'black', '7' = 'black')),
      labels = c('1','2','3','4','5','6','7'),
      labels_gp = gpar(col = 'white', fontsize = 6, fontfamily='Helvetica')
    )
  )
)

clustered_col_dend <- column_dend(pval_heatmap)
labels(clustered_col_dend)

clustered_col_dend %>% 
  lapply(labels) %>% 
  enframe() %>% 
  unnest(value) %>% 
  rename('cluster' = name, 'strain' = value) %>% 
  write_csv('Supplementary_Data_Tables/Excel_files/corr_clustering_heatmap_cluster_gene_sets.csv')

# loop to extract genes for each cluster.
 for (i in 1:length(clustered_col_dend)){
   
   
   if (i == 1) {
     clu <- t(t(row.names(mat[row_order(HM)[[i]],])))
     out <- cbind(clu, paste("cluster", i, sep=""))
     colnames(out) <- c("GeneID", "Cluster")
   } else {
   clu <- t(t(row.names(mat[row_order(HM)[[i]],])))
   clu <- cbind(clu, paste("cluster", i, sep=""))
   out <- rbind(out, clu)
   }
}


labels(clustered_col_dend[[5]])


# make legends
heatmap_legend <-
  Legend(title = 'P-value of pos. corr',
         col_fun = colorRamp2(c(0, 0.1), c(ucsf_colors$purple1, 'white')),
         at = c(0, 0.05, 0.1),
         direction = 'horizontal',
         grid_height = unit(1, 'mm'),
         grid_width = unit(0.5, "mm"),
         title_gp = gpar(fontsize = 6),
         title_position = 'topcenter',
        labels_gp = gpar(fontsize = 6))

barchart_legend <-
  Legend(title = 'log(GAP eff - GEF eff)',
         col_fun = GAPGEF_col_fn,
         at = c(-2, 0, 3.5),
         direction = 'horizontal',
         grid_height = unit(1, 'mm'),
         grid_width = unit(0.5, "mm"),
         title_gp = gpar(fontsize = 6),
         title_position = 'topcenter',
         labels_gp = gpar(fontsize = 6))

# save 
pdf('Figure4_Multispecificity/Plots/4B_Corr_Pvalue_Heatmap.pdf', width = 4.5, height = 2.7)
draw(pval_heatmap)
dev.off()

pdf('talks/Corr_Pvalue_Heatmap_for_slides.pdf', width = 6, height = 2.7)
draw(pval_heatmap)
dev.off()

pdf('Figure4_Multispecificity/Plots/4B_Corr_Pvalue_Legend.pdf', width = 1, height = 1)
draw(heatmap_legend)
dev.off()

pdf('Figure4_Multispecificity/Plots/4B_Kinetics_Legend.pdf', width = 1, height = 1)
draw(barchart_legend)
dev.off()


# also plot pcc_heatmap
pcc_corr_mat <-
  filtered_correlations %>%
  filter(grepl('GSP1', query_uniq1) & (! grepl('GSP1', query_uniq2))) %>%
  mutate('mutant' = substring(query_uniq1, first = 8)) %>%
  select(mutant, 'query' = query_uniq2, pearson) %>%
  spread(query, pearson) %>%
  column_to_rownames('mutant') %>%
  as.matrix()


# plotting heatmap
pcc_heatmap <-
  Heatmap(
    pcc_corr_mat, name = 'Pearson corr\ncoefficient', show_heatmap_legend = F,
    col = colorRamp2(c(-0.4, 0, 0.4), c(ucsf_colors$pink1, 'white', ucsf_colors$green1)),

    # rows
    row_title = 'Gsp1 point mutant', row_title_gp = gpar(fontsize = 6),
    cluster_rows = as.dendrogram(row.hc), row_split = 3, row_dend_width = unit(6, "mm"),
    row_names_side = 'left', row_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),

    # mutant cluster annotations
    left_annotation = HeatmapAnnotation(
      which = 'row',
      width = unit(2, "mm"),
      annot = anno_block(
        gp = gpar(fill = c('1' = ucsf_colors$navy1,
                           '2' = ucsf_colors$navy2,
                           '3' = ucsf_colors$navy3)),
        labels = c('1','2','3'),
        labels_gp = gpar(col = 'white', fontsize = 6, fontfamily='Helvetica')
      )
    ),

    # mutant barplot annotations
    right_annotation = HeatmapAnnotation(
      which = 'row',
      width = unit(4.5, "cm"),
      `Log ratio of\nGAP/GEF rel.\nefficiences` = anno_barplot(
        logGAPGEF, gp = gpar(limits = c(-2, 3.5), fill = GAPGEF_colors),
        axis_param = list(gp = gpar(fontsize = 6, fontfamily='Helvetica'),
                          labels_rot = 0)),
      `in vitro\nGAP relative\nefficiency\nMUT/WT` = anno_barplot(
        GAP_for_barchart, gp = gpar(limits = c(0, 5), fill = ucsf_colors$gray1),
        axis_param = list(gp = gpar(fontsize = 6, fontfamily='Helvetica'),
                          labels_rot = 0)),
      `in vitro\nGEF relative\nefficiency\nMUT/WT` = anno_barplot(
        GEF_for_barchart, gp = gpar(limits = c(0, 5), fill = ucsf_colors$gray1),
        axis_param = list(gp = gpar(fontsize = 6, fontfamily='Helvetica'),
                          labels_rot = 0)),
      annotation_name_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
      annotation_name_side = 'top',
      annotation_name_rot = 0
    ),

    # columns
    column_title = 'Saccharomyces cerevisiae genes', column_title_gp = gpar(fontsize = 6),
    cluster_columns = as.dendrogram(col.hc), column_split = 7, column_dend_height = unit(6, "mm"),
    # cluster_columns = F, column_order = col_pcoor_ordering,
    column_order = col.hc$order,
    show_column_names = F,
    column_names_gp = gpar(fontsize = 2, fontfamily='Helvetica'),
    
    # yeast gene cluster annotation
    top_annotation = HeatmapAnnotation(
      which = 'column', height = unit(2, "mm"),
      annot = anno_block(
        gp = gpar(fill = c('1' = 'black','2' = 'black', '3' = 'black', '4' = 'black',
                           '5' = 'black','6' = 'black', '7' = 'black')),
        labels = c('1','2','3','4','5','6','7'),
        labels_gp = gpar(col = 'white', fontsize = 6, fontfamily='Helvetica')
      )
    )
  )

# make legends
pcc_heatmap_legend <-
  Legend(title = 'Pearson corr. coefficient',
         col_fun = colorRamp2(c(-0.4, 0, 0.4), c(ucsf_colors$pink1, 'white', ucsf_colors$green1)),
         at = c(-0.4, 0, 0.4),
         direction = 'horizontal',
         grid_height = unit(1, 'mm'),
         grid_width = unit(0.5, "mm"),
         title_gp = gpar(fontsize = 6),
         title_position = 'topcenter',
         labels_gp = gpar(fontsize = 6))

# save
pdf('Extended_Figures/Ext_Fig8_PCC_Heatmap.pdf', width = 7.2, height = 2.7)
draw(pcc_heatmap)
dev.off()

pdf('Extended_Figures/Ext_Fig8_PCC_Legend.pdf', width = 1, height = 1)
draw(pcc_heatmap_legend)
dev.off()

pdf('Extended_Figures/Ext_Fig8_Kinetics_Legend.pdf', width = 1, height = 1)
draw(barchart_legend)
dev.off()
