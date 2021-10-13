library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(dendextend)
source('ucsf_colors.R')

# define clustering parameters
clustfn <- function(mat) {
  cormat <- cor(mat, use = "pairwise.complete.obs", method = "pearson")
  dissim <- as.dist((1 - cormat)/2)
  hc <- hclust(dissim, method = "average" )
  return(hc)
}

# read in data
residue_data <- read_tsv('Data/interface_residues_and_mutations.txt', col_types = cols())
strong_mutant_res <- c(34, 78, 79, 80, 101, 108, 112, 141, 147, 148, 157)

# convert data to matrix
mat <-
  residue_data %>% 
  mutate('yeastresnum' = factor(yeastresnum, unique(residue_data$yeastresnum))) %>% 
  filter(yeastresnum %in% strong_mutant_res) %>% 
  arrange(desc(deltarASA)) %>% 
  select(yeastresnum, partner, deltarASA) %>% 
  spread(partner, deltarASA) %>% 
  column_to_rownames('yeastresnum') %>% 
  as.matrix()

# convert to protein name format (i.e. GSP1 -> Gsp1)
colnames(mat) <- lapply(colnames(mat), function(x) paste0(substr(x, 1, 1), tolower(substr(x, 2, nchar(x)))))

# cluster rows and columns
row.hc <- clustfn(t(mat))
col.hc <- clustfn(mat)

# # use the first principal coordinate (classical multidimensional scaling) to rotate dendrogram
row_pcoor <-
  cmdscale(dist(mat), eig = T, k = 1)$point %>%
  as_tibble(rownames = 'residue') %>%
  arrange(`V1`) %>% pull(residue)
row.hc <- rotate(row.hc, row_pcoor)

col_pcoor <-
  cmdscale(dist(t(mat)), eig = T, k = 1)$point %>%
  as_tibble(rownames = 'partner') %>%
  arrange(`V1`) %>% pull(partner)
col.hc <- rotate(col.hc, col_pcoor)

### print for source file
mat %>% as_tibble() %>% 
  mutate('residue number' = rownames(mat)) %>% 
  select('residue number', everything()) %>% 
  write_tsv('Per_Figure_source_files/EDF4a.txt')
  


# plot the Heatmap
pdf('Figure1_E-MAP/Plots/1C_Interface_Heatmap.pdf', width = 2.5, height = 1.7)
draw(
  Heatmap(mat, name = 'delta rASA',
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.rect(x, y, width, height, gp = gpar(col = 'black', fill = 'transparent', lwd = 0.5))},
          col = colorRamp2(c(0, 0.7), c('white', ucsf_colors$navy1)), na_col = 'white',
          show_heatmap_legend = F,
          border = F,
          row_title = 'Gsp1 residue number', column_title = 'Gsp1 partner',
          row_title_gp = gpar(fontsize = 6), column_title_gp = gpar(fontsize = 6),
          cluster_rows = as.dendrogram(row.hc), cluster_columns = as.dendrogram(col.hc),
          row_dend_reorder = FALSE, column_dend_reorder = FALSE,
          show_column_dend = F,
          row_dend_width = unit(3, "mm"), row_names_side = 'left',
          row_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
          column_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
          column_names_rot = 45
  )
)
dev.off()

# plot the Legend
pdf('Figure1_E-MAP/Plots/1C_legend.pdf', width = 1, height = 1)
draw(Legend(title = expression(Delta*'rASA'),
            at = c(0, 0.35, 0.7),
            col_fun = colorRamp2(c(0, 0.7), c('white', ucsf_colors$navy1)),
            direction = 'horizontal',
            grid_height = unit(1, 'mm'),
            grid_width = unit(1, "mm"),
            title_gp = gpar(fontsize = 6),
            title_position = 'topcenter',
            labels_gp = gpar(fontsize = 6),
))
dev.off()


# unused code, making heatmap in ggplot
# residue_data %>% 
#   mutate('yeastresnum' = factor(yeastresnum, unique(residue_data$yeastresnum))) %>% 
#   filter(yeastresnum %in% strong_mutant_res) %>% 
#   arrange(desc(deltarASA)) %>% 
#   ggplot(aes(x = yeastresnum, y = partner, fill = deltarASA)) +
#   geom_tile(color = ucsf_colors$gray3) +
#   scale_fill_gradientn(colors = c('white', ucsf_colors$navy1), 
#                        limits = c(0, 0.7), breaks=c(0.35, 0.7), na.value = 'white') +
#   theme_bw() +
#   xlab('Gsp1 residue number') +
#   ylab("Gsp1 partner") +
#   labs(fill = expression(Delta*'rASA')) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
#         axis.text.y = element_text(size = 5),
#         axis.title = element_text(size = 5),
#         axis.ticks = element_line(size = 0.05),
#         axis.ticks.length = unit(0.05, 'cm'),
#         axis.line = element_line(size = 0.1),
#         legend.text = element_text(size = 5),
#         legend.title = element_text(size = 5),
#         legend.key.size = unit(0.2, "cm"),
#         legend.margin = margin(0, 0, 0, 0, unit = "cm"),
#         legend.position = 'right') 
# ggsave('Figure1_E-MAP/Plots/1C_Interface_Heatmap_ggplot.pdf', width = 2.5, height = 1.7)
