
library(tidyverse)
library(gplots)
library(ComplexHeatmap)
library(circlize)
library(dendextend)
source('ucsf_colors.R')
set.seed(3)

# set colors for heatmap
GI_cyan <- '#0BC3E8'
GI_yellow <- '#FDFB00'
GI_black <- '#000000'

# define clustering function
clustfn <- function(mat) {
  cormat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")
  dissim <- as.dist((1 - cormat)/2)
  hc <- hclust(dissim, method = "average" )
  return(hc)
}

# read in raw E-MAP data, process into a matrix
emap <-
  read_tsv('Data/E-MAP/gsp1_pEMAP_avg_merged_gene_names.txt') %>%
  gather(-Gene, key = strain, value = score) %>%
  separate(strain, 'library_gene', sep = ' - ', remove = TRUE) %>%
  separate(Gene, c('Gene', 'mutant'), sep = ' - ') %>%
  select(mutant, library_gene, score) %>%
  spread(library_gene, score) %>%
  filter(! mutant %in% c('NTER3XFLAG WT','CTER3XFLAG WT', 'T34N')) %>%
  column_to_rownames(var="mutant") %>%
  as.matrix()

colnames(emap) <- tolower(colnames(emap))

# cluster rows and columns
row_hc <- clustfn(t(emap))
col_hc <- clustfn(emap)

# use the first principal coordinate (classical multidimensional scaling) to flip leaves
row_pcoor <-
  cmdscale(dist(emap), eig = T, k = 1)$point %>%
  as_tibble(rownames = 'mutant') %>%
  arrange(`V1`) %>% pull(mutant)
row_hc <- rotate(row_hc, row_pcoor)

col_pcoor <-
  cmdscale(dist(t(emap)), eig = T, k = 1)$point %>%
  as_tibble(rownames = 'gene') %>%
  arrange(`V1`) %>% pull(gene)
col_hc <- rotate(col_hc, col_pcoor)

row_order <- row_hc$order
col_order <- col_hc$order


strong_mutants <- c('H141E','Y157A','D79A','D79S','T34Q', 'T34E',
                    'K101R','T34G','T34A','R108L','R108I','Q147E',
                    'H141I','G80A','H141R','R112A','R112S','R108Y',
                    'R108Q','Y148I','R108G','R108A','R78K')

# col_hc is our hierarchical clustering of array genes.
# We split this into the k=1200 most similar clusters
cluster_assignments <- cutree(col_hc, k = 1200)

# add SGD descriptions to each gene
clusters_with_descriptions <-
  data.frame('name' = names(cluster_assignments),
             'cluster_number' = unname(cluster_assignments),
             stringsAsFactors = FALSE) %>%
  left_join(read_delim('Data/E-MAP/SGD_descriptions_library_genes.txt',
                       delim='\t', col_types = cols()), by = 'name')

# write out each cluster with more than 2 members to a CSV file to examine descriptions
cluster_dir = 'Data/E-MAP/raw_data_clustering/'
dir.create(file.path(getwd(), cluster_dir))
clusters_with_descriptions %>%
  group_by(cluster_number) %>%
  filter(n() > 2) %>%
  nest(.key = 'cluster_data') %>%
  mutate(data = walk2(cluster_data,
                      paste0(cluster_dir, '/cluster_', cluster_number, '.csv'),
                      delim = ',', na = '', write_delim))

# for the main figure, only show the strong mutants
# the clustering separates these mutants out as the first branch of the dendrogram
strong_dend <- as.dendrogram(row_hc)[[1]]
order.dendrogram(strong_dend) <- seq_along(get_nodes_attr(strong_dend, 'members'))
emap_strong <- emap[labels(strong_dend),]

# get list of gene names for each cluster, gene names hand chosen based on descriptions from SGD
clusts <- list()
clusts$`mRNA\nexport` <- pull(filter(clusters_with_descriptions, cluster_number %in% c(612, 613)), name)
clusts$Dynactin <- pull(filter(clusters_with_descriptions, cluster_number %in% c(61, 262)), name)
clusts$`tRNA\nurmylation` <- c('elp2','elp6','elp3','elp4','iki1')
clusts$`Spindle assembly\nregulation` <- pull(filter(clusters_with_descriptions, cluster_number == 161), name)
clusts$SWR1 <- pull(filter(clusters_with_descriptions, cluster_number == 63), name)
clusts$HOG1 <- pull(filter(clusters_with_descriptions, cluster_number %in% c(109, 639)), name)
clusts$`mRNA\nsplicing` <- c('mud1', 'not5', 'gbp2')
clusts$Mitochondrial <- c('crd1', 'mdm38', 'rim1', 'pet123', 'mip1', 'rpo41', 'img2')
clusts$Rpd3L <- pull(filter(clusters_with_descriptions, cluster_number == 235), name)

# Make heatmap of full data, showing only strong mutants
hm <-
  Heatmap(emap_strong, name = 'GI S-score',
          col = colorRamp2(c(-4, 0, 4), c(GI_cyan, GI_black, GI_yellow)),
          bottom_annotation = columnAnnotation(
            inset = colnames(emap_strong) %in% unlist(clusts),
            show_annotation_name = F, show_legend = F,
            col = list(inset = c('TRUE' = 'black', 'FALSE' = 'white')),
            simple_anno_size = unit(2, "mm")),
          show_heatmap_legend = F,
          column_title = 'S. cerevisiae E-MAP deletion library (n = 1444)',
          row_title = 'Gsp1 point mutant',
          row_title_gp = gpar(fontsize = 7), column_title_gp = gpar(fontsize = 7),
          cluster_rows = strong_dend, cluster_columns = as.dendrogram(col_hc),
          row_dend_reorder = FALSE, column_dend_reorder = FALSE,
          row_dend_width = unit(5, "mm"), row_dend_side = "left",
          show_column_dend = F, show_column_names = F,
          row_names_side = "left", row_names_gp = gpar(fontsize = 6, fontfamily='Helvetica')
          )
### for Figure source file
emap_strong %>% as_tibble() %>% 
  mutate('Gsp1_mutant' = rownames(emap_strong)) %>% 
  select(Gsp1_mutant, everything()) %>% 
  write_csv('Per_Figure_source_files/Fig1C.csv')
emap %>% as_tibble() %>% 
  mutate('Gsp1_mutant' = rownames(emap)) %>% 
  select(Gsp1_mutant, everything()) %>% 
  write_csv('Per_Figure_source_files/EDF3.csv')

fig_slices <- c(clusts$`mRNA\nexport`, clusts$`tRNA\nurmylation`, clusts$`Spindle assembly\nregulation`)
edf_slices <- c(clusts$Dynactin, clusts$SWR1, clusts$HOG1, clusts$`mRNA\nsplicing`, clusts$Mitochondrial, clusts$Rpd3L)
emap_strong %>% as_tibble() %>% 
  mutate('Gsp1_mutant' = rownames(emap_strong)) %>% 
  select(Gsp1_mutant, fig_slices) %>% 
  write_csv('Per_Figure_source_files/Fig1C_slices.csv')
emap_strong %>% as_tibble() %>% 
  mutate('Gsp1_mutant' = rownames(emap_strong)) %>% 
  select(Gsp1_mutant, edf_slices) %>% 
  write_csv('Per_Figure_source_files/EDF4B.csv')


####
# prepare a function to plot the heatmap of a small set of a genes (a "slice")
make_heatmap_of_cluster <- function(slice, cluster_name) {
  Heatmap(slice, name = cluster_name, show_heatmap_legend = F,
          width = unit (1.6, 'mm')*ncol(slice),
          col = colorRamp2(c(-4, 0, 4), c(GI_cyan, GI_black, GI_yellow)),
          cluster_rows = strong_dend, cluster_columns = F,
          row_dend_reorder = F, column_dend_reorder = F,
          show_column_names = T, show_column_dend = F,
          show_row_dend = F, show_row_names = F, 
          column_title = cluster_name, column_title_gp = gpar(fontsize = 6),
          column_names_side = 'bottom',
          column_names_gp = gpar(fontsize = 6, fontfamily='Helvetica', fontface = 'italic')
  )
}

# add heatmaps for the small gene sets that will be shown as slices in the main figure
# main_fig_clusts <- clusts[names(clusts) %in% c('mRNA export','Dynactin','tRNA\nurmylation','Spindle assembly\nregulation')]
main_fig_clusts <- clusts[names(clusts) %in% c('mRNA\nexport','tRNA\nurmylation','Spindle assembly\nregulation')]

for (i in seq_along(main_fig_clusts)) {
  genes <- unlist(main_fig_clusts[i])
  order <- intersect(colnames(emap)[col_order], genes)
  cluster_name <- names(main_fig_clusts[i])
  slice <- emap[strong_mutants, colnames(emap) %in% genes][,order]
  hm = hm + make_heatmap_of_cluster(slice, cluster_name)
}

# plot the heatmap, with slice
# pdf('Figure1_E-MAP/Plots/1D_E-MAP.pdf', height = 2.6, width = 7.2)
# pdf('Figure1_E-MAP/Plots/1D_E-MAP_FinalFormat.pdf', height = 2.4, width = 4.9)
pdf('Figure1_E-MAP/Plots/1D_E-MAP_FinalFormat.pdf', height = 2.4, width = 7.2)
draw(hm, ht_gap = unit(0.6, 'mm'))
dev.off()

# plot the legend for the main figure heatmap
legend <- Legend(title = 'GI S-score',
                 col_fun = colorRamp2(c(-4, 0, 4), c(GI_cyan, GI_black, GI_yellow)),
                 direction = 'horizontal',
                 grid_height = unit(1, 'mm'),
                 grid_width = unit(1, "mm"),
                 title_gp = gpar(fontsize = 6),
                 title_position = 'topcenter',
                 labels_gp = gpar(fontsize = 6),
                )
pdf('Figure1_E-MAP/Plots/1D_legend.pdf', width = 1, height = 1)
draw(legend)
dev.off()


# Supplemental Figures

# Full heatmap showing weak mutants, as well as array dendrogram
full_hm <- Heatmap(emap, name = 'S-score', show_heatmap_legend = F,
                   col = colorRamp2(c(-4, 0, 4), c(GI_cyan, GI_black, GI_yellow)),
                   row_title = 'Gsp1 Point Mutant', row_title_gp = gpar(fontsize = 6),
                   cluster_rows = as.dendrogram(row_hc), cluster_columns = as.dendrogram(col_hc),
                   row_dend_reorder = FALSE, column_dend_reorder = FALSE,
                   row_dend_width = unit(6, "mm"), column_dend_height = unit(6, "mm"),
                   row_names_side = "left", row_dend_side = "left",
                   column_title = 'S. cerevisiae E-MAP deletion library (n = 1444)',
                   column_title_gp = gpar(fontsize = 6), show_column_names = F,
                   row_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'))
pdf('Extended_Figures/Ext_Fig3A_Full_E-MAP.pdf', height = 5, width = 7.2)
draw(full_hm)
dev.off()

# new function for supplemental gene set slices, which need the dendrogram
make_heatmap_of_supp_cluster <- function(slice, cluster_name) {
  if_first <- ifelse(i == 1, T, F)
  Heatmap(slice, name = cluster_name, show_heatmap_legend = F,
          width = unit (1.8, 'mm')*ncol(slice),
          col = colorRamp2(c(-4, 0, 4), c(GI_cyan, GI_black, GI_yellow)),
          cluster_rows = strong_dend, cluster_columns = F,
          row_dend_reorder = F, column_dend_reorder = F,
          show_column_names = T, show_column_dend = F,
          show_row_dend = if_first, show_row_names = if_first,
          row_title = 'Gsp1 Point Mutant', row_title_gp = gpar(fontsize = 7),
          column_title = cluster_name, column_title_gp = gpar(fontsize = 6),
          row_names_side ='left', column_names_side = 'bottom',
          row_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
          column_names_gp = gpar(fontsize = 6, fontfamily='Helvetica', fontface = 'italic')
  )
}

# plot the supplemental gene set slices
hm_supp_slices <- HeatmapList()
# supp_fig_clusts <- clusts[! names(clusts) %in% c('mRNA export','Dynactin','tRNA\nurmylation','Spindle assembly\nregulation')]
supp_fig_clusts <- clusts[! names(clusts) %in% c('mRNA\nexport','tRNA\nurmylation','Spindle assembly\nregulation')]

for (i in seq_along(supp_fig_clusts)) {
  genes <- unlist(supp_fig_clusts[i])
  order <- intersect(colnames(emap)[col_order], genes)
  cluster_name <- names(supp_fig_clusts[i])
  slice <- emap[strong_mutants, colnames(emap) %in% genes][,order]
  hm_supp_slices = hm_supp_slices + make_heatmap_of_supp_cluster(slice, cluster_name)
}

# pdf('Extended_Figures/Ext_Fig3B_Extra_Slices.pdf', height = 2.8, width = 3.5)
pdf('Extended_Figures/Ext_Fig3B_Extra_Slices_Final_Formatting.pdf', height = 2.8, width = 4)
draw(hm_supp_slices, ht_gap = unit(2, 'mm'))
dev.off()

# save a legend for the Extended_Figures too
pdf('Extended_Figures/Ext_Fig3_legend.pdf', width = 1, height = 1)
draw(legend)
dev.off()
