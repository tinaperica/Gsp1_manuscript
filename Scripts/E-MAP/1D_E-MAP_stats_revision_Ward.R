
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
  hc <- hclust(dissim, method = "ward.D2" )
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

strong_dend <- as.dendrogram(row_hc)[[1]]
order.dendrogram(strong_dend) <- seq_along(get_nodes_attr(strong_dend, 'members'))
emap_strong <- emap[labels(strong_dend),]

# Make heatmap of full data, showing only strong mutants
hm <-
  Heatmap(emap_strong, name = 'GI S-score',
          col = colorRamp2(c(-4, 0, 4), c(GI_cyan, GI_black, GI_yellow)),
          column_title = 'S. cerevisiae E-MAP deletion library (n = 1444)',
          row_title = 'Gsp1 point mutant',
          row_title_gp = gpar(fontsize = 7), column_title_gp = gpar(fontsize = 7),
          cluster_rows = strong_dend, cluster_columns = as.dendrogram(col_hc),
          row_dend_reorder = FALSE, column_dend_reorder = FALSE,
          row_dend_width = unit(6, "mm"), row_dend_side = "left",
          show_column_dend = F, show_column_names = F,
          row_names_side = "left", row_names_gp = gpar(fontsize = 6, fontfamily='Helvetica')
          )


# plot the heatmap, with slice
pdf('Revisions2/Fig_and_analysis_for_StatReviewer/1D_E-MAP_Wards.pdf', height = 2.6, width = 7.2)
draw(hm)
dev.off()
