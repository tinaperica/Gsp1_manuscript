#### plot point plot with fold change
library(tidyverse)
library(factoextra)
library(ComplexHeatmap)
library(dendextend)
library(circlize)
source('ucsf_colors.R')

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
#### for clustering reasons
data <- read_tsv('Data/APMS_data.txt') %>% 
  filter(norm == 'eqM') %>% 
  filter(! Prey_gene_name == 'GSP1') %>% 
  mutate('log2FC' = as.double(log2FC)) %>% 
  mutate('log2FC' = ifelse(log2FC < -4, -4, log2FC)) %>% 
  mutate('log2FC' = ifelse(log2FC > 4, 4, log2FC)) %>% 
  mutate('log2FC' = ifelse(is.infinite(log2FC) & log2FC < 0, -4, log2FC)) %>% 
  mutate('log2FC' = ifelse(is.infinite(log2FC) & log2FC > 0, 4, log2FC)) %>% 
  mutate('Prey_gene_name' = str_c(substr(Prey_gene_name, 1, 1), tolower(substr(Prey_gene_name, 2, nchar(Prey_gene_name))))) %>% 
  mutate('sample' = ifelse(tag == 'N', str_c('N-3xFL-', mutant), str_c(mutant, '-C-3xFL'))) %>% 
  mutate('Prey_gene_name' = case_when(is.na(Prey_gene_name) ~ PreyORF,
                                      !is.na(Prey_gene_name) ~ Prey_gene_name)) %>% 
  select(sample, tag, mutant, residue, Prey_gene_name, log2FC, adj.pvalue) %>% 
  unique()

#### Source data for EDF5 b and c
data %>% write_tsv('Per_Figure_source_files/EDF5_BC.txt')
#####

#### make a heatmap of all samples versus prey partners and cluster the samples by tag
mat <-
  data %>% 
  select(sample, Prey_gene_name, log2FC) %>%
  filter(grepl('N-3xFL', sample)) %>% 
  spread(Prey_gene_name, log2FC) %>% 
  column_to_rownames('sample') %>% 
  as.matrix()

##### ward.D2 clustering
row.dist <- get_dist(mat, method = 'euclidean') 
sum(is.infinite(row.dist)) # there are 0 Infs  so set the Infs to max or min non-inf distance
sum(is.na(row.dist))
row.hc <- hclust(row.dist, method = 'ward.D2')
ordered_samples <- row.hc$labels[row.hc$order]
col.dist <- get_dist(t(mat), method = 'euclidean')
sum(is.infinite(col.dist)) # there are 0 Infs  so set the Infs to max non-inf distance
col.hc <-  hclust(col.dist, method = 'ward.D2')
ordered_preys <- col.hc$labels[col.hc$order]

#### now make a large heatmap with all the data
pdf('Extended_Figures/Ext_Fig4_all_prey_heatmap_N.pdf', width = 4, height = 2.5)
Heatmap(
  mat, name = 'log2(MUT/WT)', show_heatmap_legend = F,
  col = colorRamp2(c(-4, 0, 4), c(ucsf_colors$pink1, 'white', ucsf_colors$blue1)),
  row_title = 'N-term tagged Gsp1 point mutant', row_title_gp = gpar(fontsize = 7),
  cluster_rows = as.dendrogram(row.hc), row_dend_width = unit(6, "mm"),
  row_names_side = 'left', row_names_gp = gpar(fontsize = 6, fontfamily = 'Helvetica'),
  column_title = 'AP-MS prey protein', column_title_gp = gpar(fontsize = 7),
  cluster_columns = as.dendrogram(col.hc), column_dend_height = unit(6, "mm"),
  column_names_gp = gpar(fontsize = 6, fontfamily = 'Helvetica'),
  show_column_names = F,
  na_col = ucsf_colors$gray3
)
dev.off()

#### make a heatmap of all samples versus prey partners and cluster the samples by tag
mat <-
  data %>% 
  select(sample, Prey_gene_name, log2FC) %>%
  filter(grepl('C-3xFL', sample)) %>% 
  spread(Prey_gene_name, log2FC) %>% 
  column_to_rownames('sample') %>% 
  as.matrix()

##### ward.D2 clustering
row.hc <- get_dist(mat, method = 'euclidean') %>% hclust(method = 'ward.D2')
ordered_samples <- row.hc$labels[row.hc$order]
col.dist <- get_dist(t(mat), method = 'euclidean')
# sum(is.na(col.dist)) # there are 11076 NA (and 39010 non-NA), so set the NAs to max distance
# col.dist[is.na(col.dist)] <- max(col.dist, na.rm = T)
col.hc <-  hclust(col.dist, method = 'ward.D2')
ordered_preys <- col.hc$labels[col.hc$order]

#### now make a large heatmap with all the data
pdf('Extended_Figures/Ext_Fig4_all_prey_heatmap_C.pdf', width = 3.2, height = 2.5)
Heatmap(
  mat, name = 'log2(MUT/WT)', show_heatmap_legend = F,
  col = colorRamp2(c(-4, 0, 4), c(ucsf_colors$pink1, 'white', ucsf_colors$blue1)),
  row_title = 'C-term tagged Gsp1 point mutant', row_title_gp = gpar(fontsize = 7),
  cluster_rows = as.dendrogram(row.hc), row_dend_width = unit(6, "mm"),
  row_names_side = 'left', row_names_gp = gpar(fontsize = 6, fontfamily = 'Helvetica'),
  column_title = 'AP-MS prey protein', column_title_gp = gpar(fontsize = 7),
  cluster_columns = as.dendrogram(col.hc), column_dend_height = unit(6, "mm"),
  column_names_gp = gpar(fontsize = 6, fontfamily = 'Helvetica'),
  show_column_names = F,
  na_col = ucsf_colors$gray3
)
dev.off()

# plot the legend for the main figure heatmap
legend <- Legend(title = 'log2(MUT/WT)',
                 col = colorRamp2(c(-4, 0, 4),c(ucsf_colors$pink1, 'white', ucsf_colors$blue1)),
                 direction = 'horizontal',
                 grid_height = unit(1, 'mm'),
                 grid_width = unit(1, "mm"),
                 title_gp = gpar(fontsize = 6),
                 title_position = 'topcenter',
                 labels_gp = gpar(fontsize = 6),
)
pdf('Extended_Figures/Ext_Fig4_legend.pdf', width = 1, height = 1)
draw(legend)
dev.off()



data <- read_tsv('Data/APMS_data.txt') %>% 
  filter(norm == 'eqM') %>% 
  filter(! Prey_gene_name == 'GSP1') %>% 
  mutate('Prey_gene_name' = str_c(substr(Prey_gene_name, 1, 1), tolower(substr(Prey_gene_name, 2, nchar(Prey_gene_name))))) %>% 
  mutate('sample' = ifelse(tag == 'N', str_c('N-3xFL-', mutant), str_c(mutant, '-C-3xFL'))) %>% 
  mutate('Prey_gene_name' = case_when(is.na(Prey_gene_name) ~ PreyORF,
                                      !is.na(Prey_gene_name) ~ Prey_gene_name)) %>% 
  select(sample, tag, mutant, residue, Prey_gene_name, log2FC, adj.pvalue) %>% 
  unique()
max_value <- data %>% filter(! abs(log2FC) == Inf) %>% pull(log2FC) %>% max(na.rm = T)
min_value <- data %>% filter(! abs(log2FC) == Inf) %>% pull(log2FC) %>% min(na.rm = T)
### set Inf and -Inf values to max and min values for calculating IQR
data <- data %>% 
  mutate('log2FC' = ifelse(log2FC == Inf, max_value, log2FC)) %>% 
  mutate('log2FC' = ifelse(log2FC == -Inf, min_value, log2FC))
interquartile_range <- data %>% 
  group_by(Prey_gene_name) %>% 
  summarize('IQR' = IQR(log2FC, na.rm = T)) %>% 
  arrange(desc(IQR))
interquartile_range_partners <- interquartile_range %>% 
  filter(Prey_gene_name %in% c('Rna1', 'Srm1', 'Yrb1', 'Kap95', 'Srp1', 'Pse1', 'Spa2', 'Pup2')) %>% 
  bind_rows(., tibble('Prey_gene_name' = c('mean', '+1 s.d.'), 
                      'IQR' = c(mean(interquartile_range$IQR), mean(interquartile_range$IQR) + sd(interquartile_range$IQR))))
interquartile_range %>% ggplot(aes(x = IQR)) + 
  geom_histogram(color = ucsf_colors$gray3, fill = ucsf_colors$gray3, alpha = 0.2) +
  theme_classic() +
  labs(color = 'Gsp1 partner as AP-MS prey protein') +
  xlab(expression('Interquartile range of log'[2]*'(abundance'[MUT]*'/abundance'[WT]*') for each prey protein')) +
  geom_segment(data = interquartile_range_partners, 
               aes(x = IQR, xend = IQR, y = 0, yend = 20, color = Prey_gene_name),
               arrow = arrow(length = unit(0.2,"cm"), ends = 'first'),
               size = 1.1, alpha = 0.8) +
  scale_color_manual(limits = c("Rna1", 'Srp1', "Kap95", "Srm1", 'Yrb1', 'Pse1', 'mean', '+1 s.d.'),
                     values=c(ucsf_colors$orange1, ucsf_colors$pink2, ucsf_colors$green1, ucsf_colors$cyan1, ucsf_colors$yellow1, ucsf_colors$green3, ucsf_colors$gray1, ucsf_colors$gray2)) +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.ticks = element_line(size = 0.05),
        axis.ticks.length = unit(0.05, 'cm'),
        legend.position = 'right',
        legend.direction = 'vertical', 
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        #legend.key.size = unit(0.3, "cm"),
        axis.line = element_line(size = 0.1))
ggsave('Extended_Figures/Ext_Fig4_IQR_of_log2FC.pdf', width = 7.2, height = 2.3) 
write_tsv(interquartile_range, 'Supplementary_Data_Tables/Supp_Table8_IQR_of_AP-MS_preys.txt')

data <- read_tsv('Data/APMS_data.txt') %>% 
  select(sample, tag, mutant, 'normalization_method' = norm, PreyORF, Prey_gene_name, log2FC, adj.pvalue, abundance) %>% 
  unique()
write_tsv(data, 'Supplementary_Data_Tables/Supp_Table7_AP-MS_data.txt')
# #### the code below makes clustered heatmaps of prays versus mutants
# ### two heatmaps: 1. is struct. partners only, second is all preys
# # core_partners
# core_partners <- data %>% pull(interface_partner) %>% unique()
# #core_partners <- append(core_partners, 'MOG1')
# ordered_mutants <- data %>% select(residue, mutant) %>% unique() %>% arrange(residue) %>% pull(mutant)
# ### preys
# preys <- data %>% pull(Prey_gene_name) %>% unique()
# 
# # make a matrix of mutants vs. partners
# mat <-
#   data %>% 
#   filter(Prey_gene_name %in% core_partners & interface_partner %in% preys & Prey_gene_name == interface_partner) %>% 
#   select(mutant, tag, Prey_gene_name, log2FC, interface_partner, deltarASA) %>% 
#   mutate('sample' = ifelse(tag == 'N', str_c('N-3xFL-', mutant), str_c(mutant, '-C-3xFL'))) %>% 
#   mutate('log2FC' = ifelse(log2FC < -4, -4, log2FC)) %>% 
#   mutate('log2FC' = ifelse(log2FC > 4, 4, log2FC)) %>% 
#   mutate('mutant' = factor(mutant, levels = ordered_mutants)) %>% 
#   arrange(mutant) %>% 
#   mutate('mutant' = paste(mutant, tag, sep = '_')) %>% 
#   rename('partner' = Prey_gene_name) %>% 
#   select(mutant, partner, log2FC) %>% 
#   spread(partner, log2FC) %>% 
#   column_to_rownames('mutant') %>% 
#   as.matrix()
# 
# 
# 
# ##### ward.D2 clustering
# 
# row.hc <- clustfn(mat, dist_method = 'euclidean', clust_method = 'ward.D2')
# col.hc <- clustfn(t(mat), dist_method = 'euclidean', clust_method = 'ward.D2')
# 
# pdf('APMS_Figure2/APMS_partners_clustered.pdf', width = 3, height = 4)
# Heatmap(
#   mat, name = 'log2FC', show_heatmap_legend = T,
#   col = colorRamp2(c(-4, 0, 4), c(ucsf_colors$pink1, 'white', ucsf_colors$blue1)),
#   
#   row_title = 'Gsp1 point mutant - tag', row_title_gp = gpar(fontsize = 8),
#   cluster_rows = as.dendrogram(row.hc), row_dend_width = unit(6, "mm"),
#   row_names_side = 'left', row_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
#   column_title = 'Core Partner', column_title_gp = gpar(fontsize = 8),
#   cluster_columns = as.dendrogram(col.hc), column_dend_height = unit(6, "mm"),
#   column_names_gp = gpar(fontsize = 6, fontfamily='Helvetica')
# )
# dev.off()
# 
# # remake plot, but for all preys
# mat <-
#   data %>% 
#   select(mutant, tag, Prey_gene_name, log2FC) %>% 
#   unique() %>% 
#   mutate('log2FC' = ifelse(log2FC < -4, -4, log2FC)) %>% 
#   mutate('log2FC' = ifelse(log2FC > 4, 4, log2FC)) %>% 
#   #mutate('mutant' = factor(mutant, levels = ordered_mutants)) %>% 
#   arrange(mutant) %>% 
#   group_by(mutant, Prey_gene_name) %>%
#   mutate(tag_avg_log2FC = mean(log2FC, na.rm=TRUE)) %>% 
#   ungroup() %>% 
#   mutate('prey' = Prey_gene_name) %>% 
#   select(mutant, prey, tag_avg_log2FC) %>% 
#   unique() %>%
#   spread(prey, tag_avg_log2FC) %>% 
#   column_to_rownames('mutant') %>% 
#   as.matrix()
# 
# row.hc <- get_dist(mat, method = 'euclidean') %>% hclust(method = 'ward.D2')
# col.dist <- get_dist(t(mat), method = 'euclidean')
# sum(is.na(col.dist)) # there is one NA, so set it to max distance
# col.dist[is.na(col.dist)] <- max(col.dist, na.rm = T)
# col.hc <-  hclust(col.dist, method = 'ward.D2')
# 
# pdf('APMS_Figure2/APMS_preys_clustered.pdf', width = 15, height = 4)
# Heatmap(
#   mat, name = 'log2FC', show_heatmap_legend = T,
#   col = colorRamp2(c(-4, 0, 4), c(ucsf_colors$pink1, 'white', ucsf_colors$blue1)),
#   
#   row_title = 'Gsp1 point mutant', row_title_gp = gpar(fontsize = 8),
#   cluster_rows = as.dendrogram(row.hc), row_dend_width = unit(6, "mm"),
#   row_names_side = 'left', row_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
#   column_title = 'Prey Gene', column_title_gp = gpar(fontsize = 8),
#   cluster_columns = as.dendrogram(col.hc), column_dend_height = unit(6, "mm"),
#   column_names_gp = gpar(fontsize = 4, fontfamily='Helvetica')
# )
# dev.off()
# 
