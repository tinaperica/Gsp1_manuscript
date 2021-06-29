
library(tidyverse)
library(gplots)
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(gridExtra)
library(grid)

source('ucsf_colors.R')
set.seed(3)



# set colors for heatmap
GI_cyan <- '#0BC3E8'
GI_yellow <- '#FDFB00'
GI_black <- '#000000'

# define clustering function
clustfn <- function(mat, clust_method='average') {
  cormat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")
  dissim <- as.dist((1 - cormat)/2)
  hc <- hclust(dissim, method = clust_method)
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

# do an initial clustering to get the order of library genes
col_hc <- clustfn(emap)
col_pcoor <-
  cmdscale(dist(t(emap)), eig = T, k = 1)$point %>%
  as_tibble(rownames = 'gene') %>%
  arrange(`V1`) %>% pull(gene)
col_hc <- rotate(col_hc, col_pcoor)
genes_order <- labels(as.dendrogram(col_hc))

make_emap_heatmap <- function(emap, name, title, clust_method='average') {
  # cluster rows and columns
  row_hc <- clustfn(t(emap), clust_method='average')
  col_hc <- clustfn(emap, clust_method='average')
  
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

  if (i == 2) {
    col.dend <- as.dendrogram(col_hc)
    col.dend <- rotate(col.dend, rev(labels(col.dend)))
    col_hc <- as.hclust(col.dend)
  }

  rownames(emap)[row_order]

  # rotate the dendrogram

  figure_order <- c('H141E','Y157A','D79A','D79S','T34Q',
                    'T34E','K101R','T34G','T34A','R108L',
                    'R108I','Q147E','H141I','G80A','H141R',
                    'R112A','R112S','R108Y','R108Q','Y148I',
                    'R108G','R108A','R78K','K129E','N84Y',
                    'A180T','T34D', 'K129F','T34Y', 'K129T',
                    'T139R','N102K','K143H','R108D','N102M',
                    'E115I','Q147L','E115A','F58A','K129I',
                    'T137G','T34S','GSP1-NAT','F58L',
                    'T139A','K169I','K132H','N105L','N105V',
                    'T34L','K143W','H141V','K143Y','R108S',
                    'K154M','N102I')

  row.dend <- as.dendrogram(row_hc)
  row.dend <- rotate(row.dend, figure_order)
  row_hc <- as.hclust(row.dend)

  col.dend <- as.dendrogram(col_hc)
  genes_left <- labels(col.dend)
  col.dend <- rotate(col.dend, genes_order[genes_order %in% genes_left])
  col_hc <- as.hclust(col.dend)

  strong_dend <- as.dendrogram(row_hc)[[1]]
  order.dendrogram(strong_dend) <- seq_along(get_nodes_attr(strong_dend, 'members'))
  emap_strong <- emap[labels(strong_dend),]

  # Make heatmap of full data, showing only strong mutants
  hm <-
    Heatmap(emap_strong,
    # Heatmap(emap, name = name,
            col = colorRamp2(c(-4, 0, 4), c(GI_cyan, GI_black, GI_yellow)),
            show_heatmap_legend = F,
            column_title = title,
            row_title = name,
            row_title_gp = gpar(fontsize = 7), column_title_gp = gpar(fontsize = 7),
            cluster_rows = strong_dend, 
            # cluster_rows = as.dendrogram(row_hc),
            cluster_columns = as.dendrogram(col_hc),
            row_dend_reorder = FALSE, column_dend_reorder = FALSE,
            row_dend_width = unit(6, "mm"), row_dend_side = "left",
            show_column_dend = F, show_column_names = F,
            row_names_side = "left", row_names_gp = gpar(fontsize = 6, fontfamily='Helvetica')
            )
    
  return(grid.grabExpr(draw(hm)))
}


# for each step i in the loop, remove i * 10% of data and cluster
heatmaps <- list()
for (i in seq(0,5)) {
  print(i)
  n <- floor((1-0.1*i) * ncol(emap))
  name <- 'Gsp1 point mutant'
  cols_subset <- sample(ncol(emap), n, replace=FALSE)
  emap_tmp <- emap[,cols_subset]
  percent <- 100-(10*i)
  title <- paste0(percent,'% of library genes, n = ', n)
  heatmaps[[i+1]] <- make_emap_heatmap(emap_tmp, name, title)
}

# plot the heatmap, with slice
pdf('Revisions2/Fig_and_analysis_for_StatReviewer/SuppFig19_1D_E-MAP_data_removal.pdf', height = 8, width = 8)
grid.draw(arrangeGrob(grobs=heatmaps[1:6], ncol=2))
dev.off()
