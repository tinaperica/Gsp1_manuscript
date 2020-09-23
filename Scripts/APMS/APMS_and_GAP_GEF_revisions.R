library(tidyverse)
library(extrafont)
library(ggrepel)
loadfonts()
source('ucsf_colors.R')
partners <- c('Srm1', 'Rna1', 'Yrb1', 'Mog1', 'Kap95', 'Kap120', 'Srp1', 'Pse1',
              "Tub2" , "Imh1",  "Vps71", "Puf6",  "Swr1",  "Vps72", "Cdc14",  "Ecm1",  "Rpl6b", "Rpl8b", "Rpp2a",
              "Slx9",  "Lsm6",  "Prp43", "Sub1",  "Caf40", "Sum1", 'Pol2', 'Spa2' )
APMS <- read_tsv('Data/APMS_data.txt') %>% 
  mutate('Prey_gene_name' = str_c(substr(Prey_gene_name, 1, 1), tolower(substr(Prey_gene_name, 2, nchar(Prey_gene_name))))) %>% 
  filter(Prey_gene_name %in% partners  & norm == 'eqM') %>% 
  select(mutant, tag, Prey_gene_name, log2FC, residue) %>%
  unique() %>% 
  mutate('sample' = ifelse(tag == 'N', str_c('N-3xFL-', mutant), str_c(mutant, '-C-3xFL'))) %>% 
  mutate('log2FC' = ifelse(log2FC < -4, -4, log2FC)) %>% 
  mutate('log2FC' = ifelse(log2FC > 4, 4, log2FC)) %>% 
  arrange(residue) %>% 
  spread(Prey_gene_name, log2FC) %>% 
  mutate('difference' = Rna1 - Srm1) %>% 
  mutate('difference2' = Srm1 - Rna1) %>% 
  mutate('ratio' = Rna1/Srm1) %>% 
  mutate('ratio' = ifelse(ratio < -5, -5, ratio)) %>% 
  mutate('ratio' = ifelse(ratio > 5, 5, ratio))
ordered_mutants <- APMS %>% pull(mutant) %>% unique()
kinetics <- read_tsv("Data/kinetics_data_relative_to_WT.txt") %>% 
  filter(measure %in% c("GEF_kcat_Km", "GAP_kcat_Km")) %>% 
  select(mutant, measure, ln_rel_to_WT) %>% 
  pivot_wider(names_from = measure, values_from = ln_rel_to_WT) %>% 
  select(mutant, 'GAP' = "GAP_kcat_Km", 'GEF' = "GEF_kcat_Km")


data <- APMS %>% 
  inner_join(., kinetics, by = 'mutant')
data_for_text <- data %>% 
  group_by(mutant) %>% 
  summarise('GAP' = max(GAP), 'GEF' = max(GEF)) 

plot <- data %>% 
  ggplot(aes(x = GEF, y = GAP)) +
  geom_point(data = tibble('GAP' = 0, 'GEF' = 0, 'mutant' = 'WT'), aes(x = GEF, y = GAP), color = ucsf_colors$navy2, size = 3) +
  geom_point(aes(shape = tag, color = difference), size = 6) +
  scale_color_gradientn(limits = c(-8, 7), breaks = c(-8, 0, 7), 
                        colors = c(ucsf_colors$orange1, 'white', ucsf_colors$cyan1)) +
  #scale_color_gradientn(limits = c(-8, 7), breaks = c(-8, 0, 7), ### for diff2
   #                     colors = c(ucsf_colors$cyan1, 'white', ucsf_colors$orange1)) +
  #scale_color_gradientn(limits = c(-5, 5), breaks = c(-5, 1, 5), ### for ratio
   #                     colors = c(ucsf_colors$orange1, 'white', ucsf_colors$cyan1)) +
  scale_shape_manual(values = c("\u25D7", "\u25D6"), guide = 'none') +
  geom_abline(slope = 1, intercept = 0, alpha = 0.1) +
  xlim(c(-6, 1.5)) + ylim(c(-6, 1.5)) +
  geom_text_repel(data = data_for_text, aes(label = mutant), size = 2.15, point.padding = 0.5, segment.size = 0.1) +
  theme_classic() +
  labs(colour = 'GAP - GEF\nAP-MS (ln(MUT/WT))') +
  xlab('Relative GEF efficiency (ln(MUT/WT))') +
  ylab('Relative GAP efficiency (ln(MUT/WT))') + 
  theme(text = element_text(size = 5),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = 0.05),
        axis.ticks.length = unit(0.05, 'cm'),
        legend.position = 'bottom', 
        legend.margin = margin(c(0, 0, 0, 0)),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(.3, "cm"),
        axis.line = element_line(size = 0.1),
        plot.margin = unit(c(0, 0, 0, 0), 'cm'))

cairo_pdf(file = 'Revisions/Extended_Figures/EDF_8/EDF_8a_kinetics_vs_AP-MS_semicircle_plot_diff.pdf', family = "Arial Unicode MS", width = 2.2, height = 2.6)
plot
dev.off()

#partner_colors <- tibble('partner' = partners, 'color' = c(ucsf_colors$cyan1, ucsf_colors$orange1, rep(ucsf_colors$green1, times = length(partners) - 2)))

### now plot individual partners
plot_partner <- function(partner) {
  data %>% 
    ggplot(aes(x = GEF, y = GAP)) +
    geom_point(data = tibble('GAP' = 0, 'GEF' = 0, 'mutant' = 'WT'), aes(x = GEF, y = GAP), color = ucsf_colors$navy2, size = 3) +
    geom_point(aes_string(shape = 'tag', color = partner), size = 6) +
    scale_color_gradientn(limits = c(-4, 4), breaks = c(-4, 0, 4), 
                          colors = c(ucsf_colors$pink2, ucsf_colors$gray3, ucsf_colors$blue1)) +
    scale_shape_manual(values = c("\u25D7", "\u25D6"), guide = 'none') +
    geom_abline(slope = 1, intercept = 0, alpha = 0.1) +
    xlim(c(-6, 1.5)) + ylim(c(-6, 1.5)) +
    geom_text_repel(data = data_for_text, aes(label = mutant), size = 2.15, point.padding = 0.5, segment.size = 0.1) +
    theme_classic() +
    labs(colour = str_c(partner, ' (log2(MUT/WT))')) +
    xlab('Relative GEF efficiency') +
    ylab('Relative GAP efficiency') + 
    theme(text = element_text(size = 7),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.title = element_text(size = 6),
          axis.ticks = element_line(size = 0.05),
          axis.ticks.length = unit(0.05, 'cm'),
          legend.position = 'bottom', 
          legend.margin = margin(c(0, 0, 0, 0)),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          legend.key.size = unit(.3, "cm"),
          axis.line = element_line(size = 0.1),
          plot.margin = unit(c(0, 0, 0, 0), 'cm'))
  
}
data <- read_tsv('Data/APMS_data.txt') %>% 
  mutate('Prey_gene_name' = str_c(substr(Prey_gene_name, 1, 1), tolower(substr(Prey_gene_name, 2, nchar(Prey_gene_name))))) %>% 
  filter(Prey_gene_name %in% partners  & norm == 'eqM') %>% 
  select(mutant, tag, Prey_gene_name, log2FC, residue, adj.pvalue) %>%
  unique() %>% 
  mutate('sample' = ifelse(tag == 'N', str_c('N-3xFL-', mutant), str_c(mutant, '-C-3xFL'))) %>% 
  mutate('log2FC' = ifelse(adj.pvalue > 0.1, 0, log2FC)) %>% 
  mutate('log2FC' = ifelse(log2FC < -4, -4, log2FC)) %>% 
  mutate('log2FC' = ifelse(log2FC > 4, 4, log2FC)) %>% 
  select(-adj.pvalue) %>% 
  arrange(residue) %>% 
  spread(Prey_gene_name, log2FC) %>% 
  inner_join(., kinetics, by = 'mutant')

### for some reason cairo_pdf doesn't work from a for loop for me...
p <- 'Srm1'
plot <- plot_partner(p)
cairo_pdf(file = 'Revisions/Extended_Figures/EDF_8/EDF_8b_kinetics_vs_AP-MS_semicircle_plot_diff.pdf', family = "Arial Unicode MS", width = 2.2, height = 2.6)
plot
dev.off()
p <- 'Rna1'
plot <- plot_partner(p)
cairo_pdf(file = str_c('Figure4_Multispecificity//APMS_and_biophysics_combined/kinetics_vs_', p, '_semicircle_plot.pdf'), family = "Arial Unicode MS", width = 2.2, height = 2.6)
plot
dev.off()
p <- 'Yrb1'
plot <- plot_partner(p)
cairo_pdf(file = str_c('Figure4_Multispecificity/APMS_and_biophysics_combined/kinetics_vs_', p, '_semicircle_plot.pdf'), family = "Arial Unicode MS", width = 2.2, height = 2.6)
plot
dev.off()
p <- 'Mog1'
plot <- plot_partner(p)
cairo_pdf(file = str_c('Figure4_Multispecificity//APMS_and_biophysics_combined/kinetics_vs_', p, '_semicircle_plot.pdf'), family = "Arial Unicode MS", width = 2.2, height = 2.6)
plot
dev.off()
p <- 'Kap95'
plot <- plot_partner(p)
cairo_pdf(file = str_c('Figure4_Multispecificity//APMS_and_biophysics_combined/kinetics_vs_', p, '_semicircle_plot.pdf'), family = "Arial Unicode MS", width = 2.2, height = 2.6)
plot
dev.off()
p <- 'Pse1'
plot <- plot_partner(p)
cairo_pdf(file = str_c('Figure4_Multispecificity/APMS_and_biophysics_combined/kinetics_vs_', p, '_semicircle_plot.pdf'), family = "Arial Unicode MS", width = 2.2, height = 2.6)
plot
dev.off()
p <- 'Kap120'
plot <- plot_partner(p)
cairo_pdf(file = str_c('Figure4_Multispecificity/APMS_and_biophysics_combined/kinetics_vs_', p, '_semicircle_plot.pdf'), family = "Arial Unicode MS", width = 2.2, height = 2.6)
plot
dev.off()
p <- 'Srp1'
plot <- plot_partner(p)
cairo_pdf(file = str_c('Figure4_Multispecificity/APMS_and_biophysics_combined/kinetics_vs_', p, '_semicircle_plot.pdf'), family = "Arial Unicode MS", width = 2.2, height = 2.6)
plot
dev.off()
p <- 'Spa2'
plot <- plot_partner(p)
cairo_pdf(file = str_c('Figure4_Multispecificity/APMS_and_biophysics_combined/kinetics_vs_', p, '_semicircle_plot.pdf'), family = "Arial Unicode MS", width = 2.2, height = 2.6)
plot
dev.off()
p <- 'Pol2'
plot <- plot_partner(p)
cairo_pdf(file = str_c('Figure4_Multispecificity/APMS_and_biophysics_combined/kinetics_vs_', p, '_semicircle_plot.pdf'), family = "Arial Unicode MS", width = 2.2, height = 2.6)
plot
dev.off()
p <- 'Vps71'
plot <- plot_partner(p)
cairo_pdf(file = str_c('Figure4_Multispecificity/APMS_and_biophysics_combined/kinetics_vs_', p, '_semicircle_plot.pdf'), family = "Arial Unicode MS", width = 2.2, height = 2.6)
plot
dev.off()

