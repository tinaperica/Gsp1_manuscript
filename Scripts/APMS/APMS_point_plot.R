#### plot point plot with fold change
library(tidyverse)
library(factoextra)
library(dendextend)
library(cowplot)
library(extrafont)
loadfonts()
source('ucsf_colors.R')

### define a clustering method to cluster preys
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


data <- read_tsv('Data/APMS_data.txt') %>% 
  filter(! Prey_gene_name == 'GSP1') %>% 
  filter(norm == 'eqM') %>% 
  mutate('Prey_gene_name' = str_c(substr(Prey_gene_name, 1, 1), tolower(substr(Prey_gene_name, 2, nchar(Prey_gene_name))))) %>% 
  mutate('interface_partner' = str_c(substr(interface_partner, 1, 1), tolower(substr(interface_partner, 2, nchar(interface_partner))))) %>% 
  mutate('Prey_gene_name' = case_when(is.na(Prey_gene_name) ~ PreyORF,
                                      !is.na(Prey_gene_name) ~ Prey_gene_name)) %>% 
  mutate('sample' = ifelse(tag == 'N', str_c('N-3xFL-', mutant), str_c(mutant, '-C-3xFL'))) %>% 
  select(sample, tag, mutant, residue, Prey_gene_name, log2FC, adj.pvalue, deltarASA, interface_partner) %>% 
  unique() %>% 
  complete(interface_partner, nesting(sample, tag, mutant, residue, Prey_gene_name, log2FC, adj.pvalue), fill = list(deltarASA = 0))
# core_partners
core_partners <- data %>% pull(interface_partner) %>% unique()
ordered_mutants <- data %>% select(residue, mutant) %>% unique() %>% arrange(residue) %>% pull(mutant)
### preys
preys <- data %>% pull(Prey_gene_name) %>% unique()


partners_data <- data %>% 
  filter(Prey_gene_name %in% core_partners & interface_partner %in% preys & Prey_gene_name == interface_partner) %>% 
  select(sample, mutant, tag, Prey_gene_name, log2FC, interface_partner, deltarASA) %>% 
  mutate('log2FC' = ifelse(log2FC < -4, -4, log2FC)) %>% 
  mutate('log2FC' = ifelse(log2FC > 4, 4, log2FC)) %>% 
  mutate('deltarASA' = ifelse(deltarASA == 0, NA, deltarASA)) %>% 
  mutate('mutant' = factor(mutant, levels = ordered_mutants)) %>% 
  arrange(mutant)
matrix_for_clustering <- partners_data %>% 
  rename('partner' = Prey_gene_name) %>% 
  select(sample, partner, log2FC) %>% 
  spread(partner, log2FC) %>% 
  column_to_rownames('sample') %>% 
  as.matrix()
col.hc <- clustfn(t(matrix_for_clustering), dist_method = 'euclidean', clust_method = 'ward.D2')
#ordered_preys <- col.hc$labels[col.hc$order]
ordered_preys <- rev(c("Yrb1", "Srm1", "Rna1", "Kap95", "Pse1", "Srp1"))
partners_data <- partners_data %>% 
  mutate('Prey_gene_name' = factor(Prey_gene_name, levels = ordered_preys)) %>% 
  arrange(mutant, Prey_gene_name)

drASA_plot <- partners_data %>% 
  ggplot(aes(y = Prey_gene_name, x = mutant, color = log2FC, size = deltarASA, shape = tag)) +
  geom_point() +
  scale_color_gradientn(limits = c(-4, 4), breaks = c(-4, 0, 4), 
                        colors = c(ucsf_colors$pink1, 'white', ucsf_colors$blue1),
                        guide = 'none') +
  scale_shape_manual(values = c("\u25D7", "\u25D6"), guide = 'none') +
  scale_size(range = c(0, 7), breaks = c(0, 0.25, 0.5, 0.75)) +
  labs(size = expression(Delta*'rASA')) +
  theme_classic() +
  xlab(element_blank()) +
  ylab('\nGsp1 partner by AP-MS') + 
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = 0.05),
        axis.ticks.length = unit(0.05, 'cm'),
        legend.position = 'bottom', 
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        axis.line = element_line(size = 0.1),
        plot.margin = unit(c(0, 0, 0, 0), 'cm'))

partners_data <- data %>% 
  filter(Prey_gene_name %in% core_partners & interface_partner %in% preys) %>% 
  mutate('adj.pvalue' = ifelse(adj.pvalue >= 0.05, 0.1, adj.pvalue)) %>% 
  select(sample, mutant, tag, Prey_gene_name, log2FC, adj.pvalue) %>% 
  unique() %>% 
  mutate('log2FC' = ifelse(log2FC < -4, -4, log2FC)) %>% 
  mutate('log2FC' = ifelse(log2FC > 4, 4, log2FC)) %>% 
  mutate('mutant' = factor(mutant, levels = ordered_mutants)) %>% 
  mutate('Prey_gene_name' = factor(Prey_gene_name, levels = ordered_preys)) %>% 
  arrange(mutant, Prey_gene_name)

pvalue_plot <- partners_data %>% 
  ggplot(aes(y = Prey_gene_name, x = mutant, color = log2FC, size = adj.pvalue, shape = tag)) +
  geom_point() +
  scale_color_gradientn(limits = c(-4, 4), breaks = c(-4, 0, 4), 
                        colors = c(ucsf_colors$pink1, 'white', ucsf_colors$blue1)) +
  scale_shape_manual(values = c("\u25D7", "\u25D6"), guide = 'none') +
  labs(color = expression('log'[2]*'(abundance'['MUT']*'/abundance'['WT']*')')) +
  scale_size("p-value\nof prey\nfold change", range = c(6.5, 3), breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1)) +
  theme_classic() +
  xlab('Gsp1 point mutant - BAIT') +
  ylab('\nGsp1 partner - PREY') + 
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = 0.05),
        axis.ticks.length = unit(0.05, 'cm'),
        legend.position = 'bottom', 
        legend.direction = 'horizontal',
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        axis.line = element_line(size = 0.1),
        plot.margin = unit(c(0, 0, 0, 0), 'cm')) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))


combined_plots <- plot_grid(drASA_plot, pvalue_plot, align = 'v', rel_heights = c(1, 1.374), nrow = 2)
cairo_pdf(filename = 'Figure2_APMS/AP-MS_semicircle_plot_partners_TEST_v2.pdf', family = "Arial Unicode MS", width = 3.8, height = 3.89)
combined_plots
dev.off()

cairo_pdf(file = 'Figure2_APMS/AP-MS_semicircle_plot_partners_for_legend.pdf', family = "Arial Unicode MS", width = 6, height = 4)
combined_plots
dev.off()


gene_sets <- read_tsv('Data/gene_sets.txt') %>% 
  select(gene_name, gene_set, Description) %>% unique()
data_gene_set <- data %>% 
  select(-interface_partner, -deltarASA) %>%
  unique() %>% 
  inner_join(., gene_sets, by = c('Prey_gene_name' = 'gene_name')) %>% 
  arrange(gene_set)
data_gene_set %>% pull(Prey_gene_name) %>% unique()




preys_C <- data %>% filter(tag == 'C') %>% pull(Prey_gene_name) %>% unique() ### 104 preys with C tagged strains
preys_N <- data %>% filter(tag == 'N') %>% pull(Prey_gene_name) %>% unique() ### 265 preys with N tagged strains
intersect_preys <- intersect(preys_C, preys_N)   ### 52 proteins with both tags, including Kap95, Mog1, Pse1, Rna1, Srm1, Yrb1 and Yrb30

intersect_preys_data <- data %>% 
  filter(Prey_gene_name %in% intersect_preys) %>% 
  select(sample, mutant, tag, Prey_gene_name, log2FC, adj.pvalue) %>% 
  mutate('adj.pvalue' = ifelse(adj.pvalue >= 0.05, 0.1, adj.pvalue)) %>% 
  unique() %>% 
  mutate('log2FC' = ifelse(log2FC < -4, -4, log2FC)) %>% 
  mutate('log2FC' = ifelse(log2FC > 4, 4, log2FC)) %>% 
  mutate('mutant' = factor(mutant, levels = ordered_mutants)) %>% 
  arrange(mutant)
matrix_for_clustering <- intersect_preys_data %>% 
  rename('partner' = Prey_gene_name) %>% 
  select(sample, partner, log2FC) %>% 
  spread(partner, log2FC) %>% 
  column_to_rownames('sample') %>% 
  as.matrix()
col.hc <- clustfn(t(matrix_for_clustering), dist_method = 'euclidean', clust_method = 'ward.D2')
ordered_preys <- col.hc$labels[col.hc$order]
intersect_preys_data <- intersect_preys_data %>% 
  mutate('Prey_gene_name' = factor(Prey_gene_name, ordered_preys)) %>% 
  arrange(mutant, Prey_gene_name)

intersect_plot <- intersect_preys_data %>% 
  ggplot(aes(x = Prey_gene_name, y = mutant, color = log2FC, size = adj.pvalue, shape = tag)) +
  geom_point() +
  scale_color_gradientn(limits = c(-4, 4), breaks = c(-4, 0, 4), 
                        colors = c(ucsf_colors$pink1, 'white', ucsf_colors$blue1)) +
  scale_shape_manual(values = c("\u25D7", "\u25D6"), guide = 'none') +
  labs(color = expression('log'[2]*'(abundance'['MUT']*'/abundance'['WT']*')')) +
  scale_size("FDR", range = c(5.2, 2.4), breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1)) +
  theme_classic() +
  ylab('Gsp1 point mutant') +
  xlab('\nIntersect of Gsp1 interaction partners identified by AP-MS of Gsp1 strains with both tags') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.ticks = element_line(size = 0.05),
        axis.ticks.length = unit(0.05, 'cm'),
        legend.position = 'bottom', 
        legend.direction = 'horizontal',
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        #legend.key.size = unit(0.3, "cm"),
        axis.line = element_line(size = 0.1),
        plot.margin = unit(c(0, 0, 0, 0), 'cm')) +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))

cairo_pdf(file = 'Figure2_APMS/AP-MS_semicircle_plot_tag_intersect_partners_20191015.pdf', family = "Arial Unicode MS", width = 7.1, height = 3.8)
intersect_plot
dev.off()

intersect_preys_data %>% filter((Prey_gene_name %in% c('Wtm1', 'Mae1') & mutant %in% c('K132H', 'R112S')) |
                                (Prey_gene_name %in% c('Pol2', 'Ura7') & mutant %in% c('T34L')))%>% 
  select(sample, Prey_gene_name, log2FC, adj.pvalue) %>% 
  unique()




