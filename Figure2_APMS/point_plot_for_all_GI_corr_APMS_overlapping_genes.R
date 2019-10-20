library(tidyverse)
library(factoextra)

source('ucsf_colors.R')
sig_gi_corr <- read_tsv('Data/all_queries.txt') %>% 
  mutate('name' = ifelse(is.na(name), query_ORF, name)) %>% 
  select(name) %>% unique()
apms <- read_tsv('Data/APMS_data.txt') %>% 
  filter(norm == 'eqM') %>% 
  mutate('Prey_gene_name' = ifelse(is.na(Prey_gene_name), PreyORF, Prey_gene_name)) %>% 
  select(tag, mutant, residue, Prey_gene_name, log2FC, adj.pvalue) %>% 
  unique() %>% 
  filter(Prey_gene_name %in% sig_gi_corr$name)
ordered_mutants <- apms %>% select(residue, mutant) %>% unique() %>% arrange(residue) %>% pull(mutant)
overlap_partners <- apms %>% pull(Prey_gene_name) %>% unique()  ### there are only 37 genes that are both in filtered_correlations and in the apms data
apms <- apms %>% 
  filter(Prey_gene_name %in% overlap_partners) %>% 
  filter(! Prey_gene_name %in% c('CFD1','SEN1', 'SUB1', 'SLX9', 'TUB2')) %>%  ## these have only Inf or NA values 
  mutate('log2FC' = ifelse(abs(log2FC) == Inf, 0, log2FC)) %>%
  select(mutant, tag, Prey_gene_name, log2FC, adj.pvalue) %>% 
  unique() %>% 
  mutate('sample' = ifelse(tag == 'N', str_c('N-3xFL-', mutant), str_c(mutant, '-C-3xFL')))

### cluster and order the preys
hc.row <- apms %>%
  select(sample, Prey_gene_name, log2FC) %>% 
  spread(Prey_gene_name, log2FC) %>% 
  select(-sample) %>% 
  as.matrix() %>% 
  t() %>% 
  get_dist(method = 'euclidean') %>% 
  hclust(method = 'single')
data <- apms %>% spread(sample, log2FC)
mat <- as.matrix(data[-1])
rownames(mat) <- data$Prey_gene_name
row_pcoor <-
    cmdscale(dist(mat), eig = T, k = 1)$point %>%
    as_tibble(rownames = 'sample') %>%
    arrange(`V1`) %>% pull(sample)
hc.row <- rotate(hc.row, row_pcoor)
preys_ordered <- hc.row$labels[hc.row$order]

partners_data <- apms %>% 
  mutate('log2FC' = ifelse(log2FC < -4, -4, log2FC)) %>% 
  mutate('log2FC' = ifelse(log2FC > 4, 4, log2FC)) %>% 
  mutate('mutant' = factor(mutant, levels = ordered_mutants)) %>% 
  mutate('Prey_gene_name' = factor(Prey_gene_name, levels = preys_ordered)) %>% 
  arrange(mutant, Prey_gene_name)
partners_N_data <- partners_data %>% 
  filter(tag == 'N') %>% 
  mutate('mutant' = factor(mutant, levels = ordered_mutants)) %>% 
  arrange(mutant)
partners_C_data <- partners_data %>% 
  filter(tag == 'C') %>% 
  mutate('mutant' = factor(mutant, levels = ordered_mutants)) %>% 
  arrange(mutant)

pvalue_plot <- partners_data %>% 
  ggplot(aes(x = Prey_gene_name, y = mutant, color = log2FC, size = adj.pvalue, shape = tag)) +
  geom_point() +
  scale_color_gradientn(limits = c(-4, 4), breaks = c(-4, 0, 4), 
                        colors = c(ucsf_colors$pink1, 'white', ucsf_colors$blue1)) +
  scale_shape_manual(values = c("\u25D7", "\u25D6"), guide = 'none') +
  labs(color = expression('log'[2]*'(abundance'['MUT']*'/abundance'['WT']*')')) +
  scale_size("FDR", range = c(5, 0), breaks = c(0, 0.001, 0.05, 0.1)) +
  theme_classic() +
  ylab('Gsp1 point mutant') +
  xlab('\nGsp1 partner by AP-MS') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
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

cairo_pdf(file = 'APMS_Figure2/apms_gi_corr_overlaping_genes_point_plot.pdf', family = "Arial Unicode MS", width = 7, height = 3.8)
pvalue_plot
dev.off()

