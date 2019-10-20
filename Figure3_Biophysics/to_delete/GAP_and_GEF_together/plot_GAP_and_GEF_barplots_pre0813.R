#### GAP and GEF assay in the same plot
library(tidyverse)
library(ggrepel)
library(cowplot)
source('ucsf_colors.R')
directory <- 'Biophysics_Figure3/GAP_and_GEF_together'
std.error.product <- function(x, sd.x, y, sd.y) {
  sqrt( (sd.x/x)^2 + (sd.y/y)^2 ) * x/y
}

GAP.data <- read_tsv("Data/RanGAP_assay/GAP_kinetics_MichaelisMenten_parameters.txt")
GEF.data <- read_tsv("Data/RanGEF_assay/GEF_kinetics_MichaelisMenten_parameters.txt")
intrinsic_hydrolysis <- read_tsv('Data/RanGAP_assay/intrinsic_hydrolysis.txt')

intersect_mutants <- intersect(GAP.data$mutant, GEF.data$mutant)

GAP.data <- GAP.data %>% 
  filter(mutant %in% intersect_mutants)
GEF.data <- GEF.data %>% 
  filter(mutant %in% intersect_mutants)
interface_mutations <- tibble('interface' = c('in GAP interface', rep('in GEF interface', 8)), 
                              'mutant' = c('K132H', 'K101R', 'R108L', 'R108I', 'R108Y', 'R018A', 'N105L', 'R108Q', 'R108G'))

GAP_interface_mutations <- c('K132H')
mut_ordered_by_kcat_Km <- GAP.data %>% select(mutant, kcat_Km) %>% arrange(kcat_Km) %>% unique() %>% pull(mutant)
GAP_barplot_colors <- c(ucsf_colors$orange1, ucsf_colors$gray2, ucsf_colors$navy2)
gap_plot <- GAP.data %>% 
  mutate("interface" = ifelse( (mutant %in% GAP_interface_mutations), "in RanGAP interface", "not in RanGAP interface" )) %>%
  mutate('interface' = ifelse( mutant == 'WT', 'WT', interface)) %>% 
  unique() %>% 
  mutate("mutant" = factor(mutant, mut_ordered_by_kcat_Km)) %>%
  ggplot(aes(mutant, kcat_Km, fill = interface)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = kcat_Km - sd, ymax = kcat_Km + sd), width = 0.5, size = 0.2, alpha = 0.75) +
  ylab(expression("k"['cat']*" / K"['m']*" [s "^-{}^{1}*mu*"M"^-{}^{1}*"]")) +
  xlab(element_blank()) +
  ggtitle('RanGAP mediated Ran:GTP to Ran:GDP + Pi hydrolysis') + 
  labs(fill = element_blank()) +
  theme_classic() +
  scale_fill_manual(values = GAP_barplot_colors) +
  geom_hline(yintercept = GAP.data$kcat_Km[GAP.data$mutant == 'WT'], linetype = "dashed", size = 0.3) +
  theme(text = element_text(size = 6), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6), 
        axis.title = element_text(size = 6.5),
        legend.text = element_text(size = 6),
        legend.key.size =  unit(0.22, 'cm'),
        axis.line = element_line(size = 0.1),
        legend.position = 'bottom',
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "gray"))

gap_kcat_plot <- GAP.data %>% 
  mutate("interface" = ifelse( (mutant %in% GAP_interface_mutations), "in RanGAP interface", "not in RanGAP interface" )) %>%
  mutate('interface' = ifelse( mutant == 'WT', 'WT', interface)) %>% 
  unique() %>% 
  mutate("mutant" = factor(mutant, mut_ordered_by_kcat_Km)) %>%
  ggplot(aes(mutant, kcat, fill = interface)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = kcat - kcat_sd, ymax = kcat + kcat_sd), width = 0.5, size = 0.2, alpha = 0.75) +
  ylab(expression("k"['cat']*" [s"^-{}^{1}*"]")) +
  xlab(element_blank()) +
  ggtitle('RanGAP mediated Ran:GTP to Ran:GDP + Pi hydrolysis') + 
  labs(fill = element_blank()) +
  theme_classic() +
  scale_fill_manual(values = GAP_barplot_colors) +
  geom_hline(yintercept = GAP.data$kcat[GAP.data$mutant == 'WT'], linetype = "dashed", size = 0.3) +
  theme(text = element_text(size = 6), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7), 
        axis.title = element_text(size = 7.5),
        legend.text = element_text(size = 6),
        legend.key.size =  unit(0.22, 'cm'),
        axis.line = element_line(size = 0.1),
        legend.position = 'none',
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "gray"))

gap_Km_plot <- GAP.data %>% 
  mutate("interface" = ifelse( (mutant %in% GAP_interface_mutations), "in RanGAP interface", "not in RanGAP interface" )) %>%
  mutate('interface' = ifelse( mutant == 'WT', 'WT', interface)) %>% 
  unique() %>% 
  mutate("mutant" = factor(mutant, mut_ordered_by_kcat_Km)) %>%
  ggplot(aes(mutant, Km, fill = interface)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = Km - Km_sd, ymax = Km + Km_sd), width = 0.5, size = 0.2, alpha = 0.75) +
  ylab(expression("K"['m']*" ["*mu*"M]")) +
  xlab('\nGsp1 point mutant') +
  ggtitle('') + 
  labs(fill = element_blank()) +
  theme_classic() +
  scale_fill_manual(values = GAP_barplot_colors) +
  geom_hline(yintercept = GAP.data$Km[GAP.data$mutant == 'WT'], linetype = "dashed", size = 0.3) +
  theme(text = element_text(size = 6), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7), 
        axis.title = element_text(size = 7.5),
        legend.text = element_text(size = 6),
        legend.key.size =  unit(0.22, 'cm'),
        axis.line = element_line(size = 0.1),
        legend.position = 'bottom',
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "gray"))



mut_ordered_by_kcat_Km <- GEF.data %>% select(mutant, kcat_Km) %>% arrange(kcat_Km) %>% unique() %>% pull(mutant)
GEF_interface_mutations <- c('K101R', 'R108L', 'R108I', 'R108Y', 'R108A', 'N105L', 'R108Q', 'R108G')
GEF_barplot_colors <- c(ucsf_colors$cyan1, ucsf_colors$gray2, ucsf_colors$navy2)
gef_plot <- GEF.data %>% 
  mutate("GEF interface" = ifelse( (mutant %in% GEF_interface_mutations), "in RanGEF interface", "not in RanGEF interface" )) %>%
  mutate('GEF interface' = ifelse( mutant == 'WT', 'WT', `GEF interface`)) %>% 
  unique() %>% 
  mutate("mutant" = factor(mutant, mut_ordered_by_kcat_Km)) %>%
  ggplot(aes(mutant, kcat_Km, fill = `GEF interface`)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = kcat_Km - kcat_Km_sd, ymax = kcat_Km + kcat_Km_sd), width = 0.5, size = 0.2, alpha = 0.75) +
  ylab(expression("k"['cat']*" / K"['m']*" [s "^-{}^{1}*mu*"M"^-{}^{1}*"]")) +
  xlab('\nGsp1 point mutant') +
  ggtitle('RanGEF mediated nucleotide exchange (Ran:GDP to Ran:mant-GTP)') + 
  labs(fill = element_blank()) +
  theme_classic() +
  scale_fill_manual(values = GEF_barplot_colors) +
  geom_hline(yintercept = GEF.data$kcat_Km[GEF.data$mutant == 'WT'], linetype = "dashed", size = 0.3) +
  theme(text = element_text(size = 6), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6), 
        axis.title = element_text(size = 6.5),
        legend.text = element_text(size = 6),
        legend.key.size =  unit(0.22, 'cm'),
        legend.position = 'bottom',
        axis.line = element_line(size = 0.1),
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "gray"))

gef_kcat_plot <- GEF.data %>% 
  mutate("GEF interface" = ifelse( (mutant %in% GEF_interface_mutations), "in RanGEF interface", "not in RanGEF interface" )) %>%
  mutate('GEF interface' = ifelse( mutant == 'WT', 'WT', `GEF interface`)) %>% 
  unique() %>% 
  mutate("mutant" = factor(mutant, mut_ordered_by_kcat_Km)) %>%
  ggplot(aes(mutant, kcat, fill = `GEF interface`)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = kcat - kcat_sd, ymax = kcat + kcat_sd), width = 0.5, size = 0.2, alpha = 0.75) +
  ylab(expression("k"['cat']*" [s"^-{}^{1}*"]")) +
  xlab('') +
  ggtitle('RanGEF mediated nucleotide exchange (Ran:GDP to Ran:mant-GTP)') + 
  labs(fill = element_blank()) +
  theme_classic() +
  scale_fill_manual(values = GEF_barplot_colors) +
  geom_hline(yintercept = GEF.data$kcat[GEF.data$mutant == 'WT'], linetype = "dashed", size = 0.3) +
  theme(text = element_text(size = 6), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7), 
        axis.title = element_text(size = 7.5),
        legend.text = element_text(size = 6),
        legend.key.size =  unit(0.22, 'cm'),
        legend.position = 'none',
        axis.line = element_line(size = 0.1),
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "gray"))
gef_Km_plot <- GEF.data %>% 
  mutate("GEF interface" = ifelse( (mutant %in% GEF_interface_mutations), "in RanGEF interface", "not in RanGEF interface" )) %>%
  mutate('GEF interface' = ifelse( mutant == 'WT', 'WT', `GEF interface`)) %>% 
  unique() %>% 
  mutate("mutant" = factor(mutant, mut_ordered_by_kcat_Km)) %>%
  ggplot(aes(mutant, Km, fill = `GEF interface`)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = Km - Km_sd, ymax = Km + Km_sd), width = 0.5, size = 0.2, alpha = 0.75) +
  ylab(expression("K"['m']*" ["*mu*"M]")) +
  xlab('\nGsp1 point mutant') +
  ggtitle('') + 
  labs(fill = element_blank()) +
  theme_classic() +
  scale_fill_manual(values = GEF_barplot_colors) +
  geom_hline(yintercept = GEF.data$Km[GEF.data$mutant == 'WT'], linetype = "dashed", size = 0.3) +
  theme(text = element_text(size = 6), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7), 
        axis.title = element_text(size = 7.5),
        legend.text = element_text(size = 6),
        legend.key.size =  unit(0.22, 'cm'),
        legend.position = 'bottom',
        axis.line = element_line(size = 0.1),
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "gray"))

gef_Km_plot_zoom <- GEF.data %>% 
  filter(! mutant %in% c('K101R', 'R108L', 'R108I', 'R108Y')) %>% 
  mutate("GEF interface" = ifelse( (mutant %in% GEF_interface_mutations), "in RanGEF interface", "not in RanGEF interface" )) %>%
  mutate('GEF interface' = ifelse( mutant == 'WT', 'WT', `GEF interface`)) %>% 
  unique() %>% 
  mutate("mutant" = factor(mutant, mut_ordered_by_kcat_Km)) %>%
  ggplot(aes(mutant, Km, fill = `GEF interface`)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = Km - Km_sd, ymax = Km + Km_sd), width = 0.5, size = 0.2, alpha = 0.75) +
  ylab(expression("K"['m']*" ["*mu*"M]")) +
  xlab('') +
  ggtitle('') + 
  labs(fill = element_blank()) +
  theme_classic() +
  scale_fill_manual(values = GEF_barplot_colors) +
  geom_hline(yintercept = GEF.data$Km[GEF.data$mutant == 'WT'], linetype = "dashed", size = 0.3) +
  theme(text = element_text(size = 6), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6), 
        axis.title = element_text(size = 7.5),
        legend.text = element_text(size = 6),
        legend.key.size =  unit(0.22, 'cm'),
        legend.position = 'none',
        axis.line = element_line(size = 0.1),
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "gray"))



### main figure - kcat/Km
plot_grid(gap_plot, gef_plot, align = 'v', nrow = 2)
ggsave(file.path(directory, 'gap_and_gef_kcat_Km_barplots.pdf'), width = 3.2, height = 5)




### supplementary figure - separate kcat and Km
plot_grid(gap_kcat_plot, gef_kcat_plot, gap_Km_plot, gef_Km_plot, align = c('v', 'h'), nrow = 2) +
  draw_plot(gef_Km_plot_zoom, x = 0.65, y = 0.22, width = 0.34, height = 0.32)
ggsave(file.path(directory, 'gap_and_gef_kcat_and_Km_supp_barplots.pdf'), width = 7.2, height = 5.3)

##### intrinsic hydrolysis barplot for supplementary data
intrinsic_hydrolysis <- intrinsic_hydrolysis %>% 
  filter(mutant %in% intersect_mutants)
mut_ordered_by_int <- intrinsic_hydrolysis %>% 
  arrange(mean_rel_rate) %>% pull(mutant)
intrinsic_hydrolysis %>% 
  mutate('residue' = ifelse( mutant == 'WT', 'WT', 'mutant')) %>% 
  mutate("mutant" = factor(mutant, mut_ordered_by_int)) %>%
  ggplot(aes(mutant, mean_rel_rate, fill = residue)) + 
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = mean_rel_rate - sd_rel_rate, ymax = mean_rel_rate + sd_rel_rate), width = 0.5, size = 0.2, alpha = 0.75) +
  ylab(expression("GTP hydrolysis rate [s"^-{}^{1}*"]")) +
  ggtitle('Intrinsic GTP hydrolysis rate') +
  xlab('\nGsp1 point mutant') +
  theme_classic() +
  scale_fill_manual(values = c(ucsf_colors$gray2, ucsf_colors$navy2)) +
  geom_hline(yintercept = intrinsic_hydrolysis$mean_rel_rate[intrinsic_hydrolysis$mutant == 'WT'], linetype = "dashed", size = 0.3) +
  theme(text = element_text(size = 6), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7), 
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.position = 'none',
        axis.line = element_line(size = 0.1),
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "gray"))
ggsave(file.path(directory, 'intrinsic_hydrolysis_barplot.pdf'), height = 3, width = 4)



#### relative GAP and GEF data on the same barplot
WT_GAP <- GAP.data %>% 
  filter(mutant == "WT") %>% 
  select(mutant, 'WT_kcat_Km' = kcat_Km, 'WT_sd' = sd)
rel_GAP.data <- GAP.data %>% 
  mutate('WT_kcat_Km' = WT_GAP$WT_kcat_Km, 'WT_sd' = WT_GAP$WT_sd) %>% 
  mutate('rel_GAP_kcat_Km' = kcat_Km/WT_kcat_Km, 
         'rel_GAP_sd' = std.error.product(x = kcat_Km, sd.x = sd, y = WT_kcat_Km, sd.y = WT_sd)) %>% 
  select(mutant, rel_GAP_kcat_Km, rel_GAP_sd)
WT_GEF <- GEF.data %>% 
  filter(mutant == "WT") %>% 
  select(mutant, 'WT_kcat_Km' = kcat_Km, 'WT_sd' = kcat_Km_sd)
rel_GEF.data <- GEF.data %>% 
  mutate('WT_kcat_Km' = WT_GEF$WT_kcat_Km, 'WT_sd' = WT_GEF$WT_sd) %>% 
  mutate('rel_GEF_kcat_Km' = kcat_Km/WT_kcat_Km, 
         'rel_GEF_sd' = std.error.product(x = kcat_Km, sd.x = kcat_Km_sd, y = WT_kcat_Km, sd.y = WT_sd)) %>% 
  select(mutant, rel_GEF_kcat_Km, rel_GEF_sd)
barplot_colors <- c(ucsf_colors$orange1, ucsf_colors$cyan1)
data <- inner_join(rel_GAP.data, rel_GEF.data, by = 'mutant') %>% 
  mutate('GAP_GEF_ratio' = rel_GAP_kcat_Km/rel_GEF_kcat_Km) %>% 
  gather(key = 'measure', value = 'value', -mutant, -rel_GAP_sd, -rel_GEF_sd, -GAP_GEF_ratio) %>% 
  mutate('sd' = ifelse(measure == 'rel_GAP_kcat_Km', rel_GAP_sd, rel_GEF_sd)) %>% 
  arrange(GAP_GEF_ratio)
order_of_mutants <- data %>% pull(mutant) %>% unique()
data <- data %>% 
  mutate('mutant' = factor(mutant, order_of_mutants)) %>% 
  select(mutant, measure, value, sd)

data %>% 
  ggplot(aes(x = mutant, y = value, fill = measure)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width = 0.3, size = 0.2, alpha = 0.75, position = position_dodge(0.9)) +
  theme_classic() +
  labs(fill = element_blank()) +
  xlab('Gsp1 point mutant') +
  ylab(label = expression("k"['cat']*"/K"['m (MUTANT)']*" / k"['cat']*"/K"['m (WILD TYPE)'])) +
  scale_fill_manual(values = barplot_colors, labels = c('GAP mediated GTP hydrolysis', 'GEF mediated nucleotide exchange')) +
  theme(text = element_text(size = 16), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 18), 
        legend.text = element_text(size = 15),
        axis.line = element_line(size = 0.1),
        legend.position = c(0.25, 0.82),
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "gray"))

ggsave(filename = file.path(directory, "combined_barplot_kcat_over_Km.pdf"), height = 5, width = 9)

