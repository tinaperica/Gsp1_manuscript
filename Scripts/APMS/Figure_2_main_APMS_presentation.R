### plot the AP-MS for main text Fig. 2
library(tidyverse)
library(ggpubr)
library(ggforce)
source('ucsf_colors.R')
interfaces <- read_tsv('Data/Gsp1_interfaces_SASA_and_conservation.txt') %>% 
  filter(protein == 'GSP1') %>% 
  mutate('interface_partner' = str_c(substr(partner, 1, 1), tolower(substr(partner, 2, nchar(partner))))) %>% 
  select(interface_partner, interface, deltarASA, deltaASA, rASAc, 'residue' = yeast_num)

apms_data <- read_tsv('Data/APMS_data.txt') %>% 
  filter(! Prey_gene_name == 'GSP1') %>% 
  filter(norm == 'eqM') %>% 
  mutate('Prey_gene_name' = str_c(substr(Prey_gene_name, 1, 1), tolower(substr(Prey_gene_name, 2, nchar(Prey_gene_name))))) %>% 
  mutate('interface_partner' = str_c(substr(interface_partner, 1, 1), tolower(substr(interface_partner, 2, nchar(interface_partner))))) %>% 
  mutate('Prey_gene_name' = case_when(is.na(Prey_gene_name) ~ PreyORF,
                                      !is.na(Prey_gene_name) ~ Prey_gene_name)) %>% 
  mutate('sample' = ifelse(tag == 'N', str_c('N-3xFL-', mutant), str_c(mutant, '-C-3xFL'))) %>% 
  select(sample, tag, mutant, residue, Prey_gene_name, log2FC, adj.pvalue, interface_partner) %>% 
  unique() %>% 
  complete(interface_partner, nesting(sample, tag, mutant, residue, Prey_gene_name, log2FC, adj.pvalue), fill = list(deltarASA = 0))

data <- apms_data %>% 
  inner_join(interfaces, by = c('interface_partner', 'residue')) %>% 
  arrange(sample, interface_partner, residue, Prey_gene_name) %>% 
  complete(interface_partner, nesting(sample, tag, mutant, residue, Prey_gene_name, log2FC, adj.pvalue), fill = list(deltarASA = 0)) %>% 
  mutate('interface' = ifelse((interface == 'rim' | interface == 'support'), 'rim/support', interface)) %>% 
  mutate('interface' = ifelse(is.na(interface), 'not_interface', interface)) 

### bin the adjusted p-value for ease of plotting - makes it more straightforward to control point size
data <- data %>% 
  mutate(pvalue_point_size = case_when(adj.pvalue <= 0.05 ~ 1,
                                   adj.pvalue > 0.05 ~ 0.1))

# core_partners
core_partners <- apms_data %>% pull(interface_partner) %>% unique()
ordered_mutants <- apms_data %>% select(residue, mutant) %>% unique() %>% arrange(residue) %>% pull(mutant)
### preys
preys <- apms_data %>% pull(Prey_gene_name) %>% unique()

data_with_interfaces <- data %>% 
  filter(interface_partner == Prey_gene_name) %>% 
  mutate('interface' = ifelse(interface == 'core', 'interface', 'not_interface'))
# yesno_interface_temp <- data_with_interfaces %>% 
#   mutate('bin_int' = ifelse(interface == 'core', 'interface', 'not_interface')) %>% 
#   filter(bin_int == 'interface') %>% 
#   select(-interface) %>% 
#   rename('interface' = bin_int)
# data_with_interfaces <- data_with_interfaces %>% bind_rows(., yesno_interface_temp)
my_comparisons <- list( c('interface', 'not_interface') )

data_with_interfaces %>%
  mutate('interface' = factor(interface, c('interface', 'not_interface'))) %>% 
  ggplot(aes(x = interface, y = log2FC, color = interface, size = pvalue_point_size)) +
  geom_sina(alpha = 1, maxwidth = 0.5, seed = 2) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = ucsf_colors$gray1, size = 0.1) +
  scale_x_discrete(breaks=c(T, F),
                   # labels = c('mutant in\npartner interface', 'mutant not in\npartner interface')
                   ) +
  scale_size(name = 'P-value of prey fold change', range = c(0.1, 1),
             breaks = c(0, 0.1, 1),
             labels = c('0', '<=0.05', '>0.05')) +
  scale_color_manual(name = 'category', guide = 'none',
                     values = c(ucsf_colors$gray1, ucsf_colors$gray2)) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar",
               color = ucsf_colors$pink1, size = 0.5, fatten = 0, width = 0.7) +
  xlab('') +
  ylab(expression('log'[2]*'FC')) +
  stat_compare_means(comparisons = my_comparisons,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                        symbols = c("****", "***", "**", "*", "ns"))) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica", size = 6),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.ticks = element_line(size = 0.05),
    axis.ticks.length = unit(0.05, 'cm'),
    legend.position = 'right',
    axis.line = element_line(size = 0.1)
  )
ggsave('Revisions/Main Figures/Figure2/APMS_sinaplot_simple_mean.pdf', height = 1.9, width = 3.3)
dev.off()
t.test(x = data_with_interfaces$log2FC[data_with_interfaces$interface == 'interface'],
    y = data_with_interfaces$log2FC[data_with_interfaces$interface == 'not_interface'])
# t = -4.6415, df = 70.442, p-value = 1.557e-05
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -2.527636 -1.008382
# sample estimates:
#   mean of x  mean of y 
# -1.0333136  0.7346954 

# t.test(x = data_with_interfaces$log2FC[data_with_interfaces$interface == 'core'],
#        y = data_with_interfaces$log2FC[data_with_interfaces$interface == 'not_interface'])
# t = -5.1099, df = 84.155, p-value = 1.984e-06
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -2.866766 -1.260580
# sample estimates:
#   mean of x mean of y 
# -1.033314  1.030359 



gef_int_resi <- interfaces %>% filter(interface_partner == 'Srm1' & interface == 'core') %>% pull(residue) %>% unique()
gap_int_resi <- interfaces %>% filter(interface_partner == 'Rna1' & interface == 'core') %>% pull(residue) %>% unique()

temp_gef <- data %>% filter(Prey_gene_name == 'Srm1') %>% 
  mutate('interface' = ifelse(residue %in% gef_int_resi, 'GEF_interface', 'not_GEF_int')) %>% 
  select(sample, residue, Prey_gene_name, log2FC, pvalue_point_size, interface) %>% 
  unique()
temp_gap <- data %>% filter(Prey_gene_name == 'Rna1') %>% 
  mutate('interface' = ifelse(residue %in% gap_int_resi, 'GAP_interface', 'not_GAP_int')) %>% 
  select(sample, residue, Prey_gene_name, log2FC, pvalue_point_size, interface) %>% 
  unique()
gap_gef <- bind_rows(temp_gef, temp_gap) %>% 
  mutate('interface' = ifelse(residue == 34, 'T34 mutation', interface)) %>% 
  mutate('interface' = factor(interface, c('GAP_interface', 'not_GAP_int', 'GEF_interface', 'not_GEF_int', 'T34 mutation'))) %>% 
  arrange(sample, Prey_gene_name, interface)
gap_gef %>% 
  ggplot(aes(x = Prey_gene_name, y = log2FC, color = interface, size = pvalue_point_size)) +
  geom_sina(alpha = 1, maxwidth = 0.5, seed = 2) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = ucsf_colors$gray1, size = 0.1) +
  scale_size(name = 'P-value of prey fold change', range = c(0.1, 1),
             breaks = c(0, 0.1, 1),
             labels = c('0', '<=0.05', '>0.05')) +
  scale_color_manual(name = 'main regulators GAP and GEF',
                     values = c(ucsf_colors$orange1, ucsf_colors$gray1, ucsf_colors$cyan1, ucsf_colors$gray1, ucsf_colors$pink1)) +
  xlab('') +
  ylab(expression('log'[2]*'FC')) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica", size = 6),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.ticks = element_line(size = 0.05),
    axis.ticks.length = unit(0.05, 'cm'),
    legend.position = 'right',
    axis.line = element_line(size = 0.1)
  )
ggsave('Revisions/Main Figures/Figure2/APMS_GAP_GEF_sinaplot.pdf', height = 2, width = 3)

# load('Data/spitzemapko_correlations_and_bonferroni_fdr_all.RData')
# correlations <- correlations %>%
#   filter(grepl('GSP1', query_uniq1), !grepl('GSP1', query_uniq2)) %>%
#   separate(query_uniq1, sep = ' - ', into = c('GSP1','query_uniq1')) %>% 
#   select('mutant' = query_uniq1, 'strain' = query_uniq2, pearson, fdr)
# name2ORF <- read_tsv('Data/spitzemap_name2ORF_index.txt') 
# SGD_descriptions <- read_tsv('Data/E-MAP/SGD_descriptions_all.txt') %>% 
#   select('partner' = name, ORF) %>% 
#   mutate('partner' = str_c(substr(partner, 1, 1), tolower(substr(partner, 2, nchar(partner)))))
# correlations <- correlations %>% 
#   inner_join(., name2ORF, by = c('strain' = 'name')) %>% 
#   inner_join(., SGD_descriptions, by = 'ORF') %>% 
#   mutate(binned_pvalue = case_when(fdr < 0.05 ~ 0.01, 
#                                    fdr <= 0.1 ~ 0.10, 
#                                    fdr > 0.1 ~ 1.00)) %>% 
#   select(mutant, 'emap_strain' = strain, 'emap_partner' = partner, 'emap_corr' = pearson, 
#          'emap_adj.pvalue' = fdr, 'emap_binned_pvalue' = binned_pvalue)
# 
# merged_apms_emap <- data %>% 
#   select(sample, mutant, Prey_gene_name, log2FC, adj.pvalue, binned_pvalue, interface) %>% 
#   inner_join(., correlations, by = c('mutant', 'Prey_gene_name' = 'emap_partner')) %>% 
#   filter(! (log2FC < -10 | log2FC > 15)) %>% 
#   mutate('point_size' = ifelse( (binned_pvalue == 0.01 & emap_binned_pvalue == 0.01), 0.1, 1 ))
# core_partners <- merged_apms_emap %>% 
#   filter(Prey_gene_name %in% c('Srm1', 'Rna1', 'Yrb1', 'Srp1', 'Kap95', 'Pse1')) %>% 
#   mutate('Prey_gene_name' = factor(Prey_gene_name, c('Srm1', 'Rna1', 'Yrb1', 'Srp1', 'Kap95', 'Pse1'))) %>% 
#   mutate('emap_strain' = factor(emap_strain, c("srm1-g282s", "srm1-ts", "rna1-1", "rna1-s116f", "yrb1-51", "srp1-5001", "kap95-e126k")))
# merged_apms_emap %>% 
#   ggplot(aes(x = log2FC, y = emap_corr, size = point_size)) + 
#   geom_point(color = ucsf_colors$gray3) +
#   geom_point(data = core_partners, aes(x = log2FC, y = emap_corr, color = emap_strain), size = 1) +
#   scale_color_manual(name = 'FDR of emap corr',
#                      values = c(ucsf_colors$cyan1, ucsf_colors$cyan3, 
#                                 ucsf_colors$orange1, ucsf_colors$orange3,
#                                 ucsf_colors$yellow2, 
#                                 ucsf_colors$pink1, ucsf_colors$green1)) +
#   scale_size(name = 'p-value of prey fold change', range = c(1, 0.1), breaks = c(0.01, 1)) +
#   ylab('Pearson correlation with partner profile') +
#   xlab(expression('log'[2]*'FC')) +
#   theme_classic() +
#   theme(
#     text = element_text(family = "Helvetica", size = 6),
#     axis.title = element_text(size = 6),
#     axis.text = element_text(size = 6),
#     axis.ticks = element_line(size = 0.05),
#     axis.ticks.length = unit(0.05, 'cm'),
#     legend.position = 'right',
#     axis.line = element_line(size = 0.1)
#   )
# ggsave('Revisions/Rebuttal_Letter_only_figures/APMS_vs_EMAP_scatterplot.pdf', height = 2.8, width = 4)
# 
# 
# 
