
##### load libraries
library(tidyverse)
library(ggforce)
library(ggrepel)
library(cowplot)
library(grid)
library(gridExtra)
library(gtable)
library(Hmisc)
source('ucsf_colors.R')

##### load datafiles
load('Data/spitzemapko_correlations_and_bonferroni_fdr_all.RData')
load('Data/filtered_v6_correlations.RData')
name2ORF <- read_tsv('Data/spitzemap_name2ORF_index.txt') # in cjm_corr

##### Make index files for Gsp1 mutants and partners
# Define strong mutants by making it past correlation filtering
mutant_index <- name2ORF %>%
  filter(grepl('GSP1', name)) %>%
  mutate('is_strong' = case_when(name %in% filtered_correlations$query_uniq1 ~ T, T ~ F)) %>%
  rename('strain' = name) %>%
  separate(strain, sep = ' - ', into = c('GSP1','mutant')) %>%
  select(-GSP1, -ORF) %>%
  filter(! mutant %in% c('GSP1-NAT','NTER3XFLAG WT', 'CTER3XFLAG WT')) %>%
  mutate('yeastresnum' = as.numeric(str_sub(mutant, 2, nchar(mutant)-1)))
strong_mutants <- select(mutant_index, mutant, is_strong)

partner_index <- data.frame(
  name = c('MSN5','SRP1','LOS1','YRB1','YRB2','KAP95','RNA1',
           'SRM1','MTR10','PSE1','NTF2','CRM1','CSE1','KAP104'),
  ORF = c('YDR335W','YNL189W','YKL205W','YDR002W','YIL063C','YLR347C','YMR235C',
          'YGL097W','YOR160W','YMR308C','YER009W','YGR218W','YGL238W','YBR017C'))
partner_strains <- filter(name2ORF, str_detect(ORF, paste(partner_index$ORF, collapse = '|')))
partner_index <-
  left_join(partner_index, partner_strains, by = 'ORF') %>%
  rename('name' = 'name.x', 'strain' = 'name.y') %>%
  filter(strain %in% correlations$query_uniq1)

core_res_table <-
  read_tsv('Data/SASA_interfaces.txt', col_types = cols()) %>%
  left_join(partner_index, by = c('partner' = 'name')) %>%
  left_join(mutant_index, by = 'yeastresnum') %>%
  filter(mutant %in% mutant_index$mutant) %>%
  mutate('is_core' = case_when(interface == 'core' ~ T, interface != 'core' ~ F)) %>%
  mutate('strain' = case_when(is.na(strain) ~ partner, T ~ strain)) %>%
  filter(is_core) %>%
  select(strain, mutant, is_core)

# clean correlations dataset and add info on whether mutation is strong
corr_for_sina <-
  correlations %>%
  filter(grepl('GSP1', query_uniq1), !grepl('GSP1', query_uniq2)) %>%
  separate(query_uniq1, sep = ' - ', into = c('GSP1','query_uniq1')) %>%
  filter(query_uniq1 %in% mutant_index$mutant) %>% 
  select(query_uniq1, query_uniq2, pearson, greater_fdr) %>%
  rename('mutant' = query_uniq1, 'strain' = query_uniq2) %>%
  left_join(strong_mutants) %>% 
  left_join(core_res_table) %>% 
  mutate(is_core = case_when(is_core ~ is_core, is.na(is_core) ~ F),
         partner = case_when(strain %in% partner_index$strain ~ T, T ~ F))

corr_for_sina %>% 
  filter(is_strong, partner) %>% 
  ggplot(aes(x = is_core, y = pearson, color = is_core, size = greater_fdr)) +
  geom_sina(alpha = 1, maxwidth = 0.8, seed = 2) +
  scale_x_discrete(breaks=c(T, F),
                   labels = c(paste('Strong mutant','in partner interface', sep = '\n'),
                              paste('Strong mutant','not in partner interface', sep = '\n'))) +
  scale_size(name = 'P-value', range = c(0.5, 0.01), breaks = c(1.0, 0.1, 0.01)) +
  scale_color_manual(name = 'category', guide = 'none',
                     values = c(ucsf_colors$gray3, ucsf_colors$gray1)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult=1), color = ucsf_colors$pink1,
               geom = "errorbar", width = 0.1, size = 0.75) +
  stat_summary(fun.y = mean, geom = "point", color = ucsf_colors$pink1, size = 2) +
  ylim(c(-0.15, 0.5)) + xlab('') + ylab('Pearson correlation\nwith partner profile') +
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
ggsave('Figure1_E-MAP/Plots/1F_Sinaplots.pdf', height = 1.9, width = 3.1)
dev.off()


corr_for_sina %>% 
  filter(is_strong, partner) %>% 
  group_by(is_core) %>% 
  summarise(n())

# Supplementary figure, showing distributions of partner/non partner for strong
data <-
  corr_for_sina %>% 
  mutate(group = case_when(is_strong & !is_core & partner ~ '1',
                           is_strong & is_core & partner ~ '2',
                           !is_strong & !is_core & partner ~ '3',
                           !is_strong & is_core & partner ~ '4',
                           !partner ~ '5'))

ggplot(data = data, aes(x = group, y = pearson, color = group)) +
  geom_sina(data = filter(data, group != '5'), aes(size = greater_fdr),
            alpha = 1, maxwidth = 0.8, seed = 2) +
  geom_violin(data = filter(data, group == '5'), fill = ucsf_colors$cyan1) +
  scale_x_discrete(breaks=c('1', '2', '3', '4', '5'),
                   labels = c(paste('Strong mutant','not in partner', 'interface', sep = '\n'),
                              paste('Strong mutant','in partner', 'interface', sep = '\n'),
                              paste('Weak mutant','not in partner', 'interface', sep = '\n'),
                              paste('Weak mutant','in partner', 'interface', sep = '\n'),
                              paste('All other','correlations', sep = '\n'))) +
  scale_size(name = 'P-value', range = c(0.5, 0.01), breaks = c(1.0, 0.1, 0.01)) +
  scale_color_manual(name = 'category', guide = 'none',
                     values = c(ucsf_colors$gray3,
                                ucsf_colors$gray1,
                                ucsf_colors$purple1,
                                ucsf_colors$purple3,
                                ucsf_colors$cyan1)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult=1), color = ucsf_colors$pink1,
               geom = "errorbar", width = 0.1, size = 0.75) +
  stat_summary(fun.y = mean, geom = "point", color = ucsf_colors$pink1, size = 2) +
  ylim(c(-0.15, 0.5)) + xlab('') + ylab('Pearson correlation') +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica", size = 6),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.ticks = element_line(size = 0.05),
    axis.ticks.length = unit(0.05, 'cm'),
    legend.position = 'top',
    legend.text = element_text(size = 6),
    axis.line = element_line(size = 0.1)
  )
ggsave('Extended_Figures/Ext_Fig3_Sinaplots_strong_weak_all.pdf', height = 2, width = 3.7)
dev.off()

data2 <-
  corr_for_sina %>% 
  mutate(group = case_when(is_strong & partner ~ '1',
                           is_strong & !partner ~ '2')) %>% 
  filter(group %in% c('1', '2'))

ggplot(data = data2, aes(x = group, y = pearson, color = group)) +
  geom_sina(data = filter(data2, group == '1'), aes(size = greater_fdr),
            alpha = 1, maxwidth = 0.8, seed = 2) +
  geom_violin(data = filter(data2, group == '2'), fill = 'black') +
  scale_x_discrete(breaks=c('1', '2'),
  labels = c(paste('Strong mutants','and partners', sep = '\n'),
             paste('Strong mutants','and non-partners', sep = '\n'))) +
  scale_size(name = 'P-value', range = c(0.5, 0.01), breaks = c(1.0, 0.1, 0.01)) +
  scale_color_manual(name = 'category', guide = 'none',
                     values = c('black', 'black')) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult=1), color = ucsf_colors$pink1,
               geom = "errorbar", width = 0.1, size = 0.75) +
  stat_summary(fun.y = mean, geom = "point", color = ucsf_colors$pink1, size = 2) +
  ylim(c(-0.15, 0.5)) + xlab('') + ylab('Pearson correlation') +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica", size = 6),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.ticks = element_line(size = 0.05),
    axis.ticks.length = unit(0.05, 'cm'),
    legend.position = 'top',
    legend.text = element_text(size = 6),
    axis.line = element_line(size = 0.1)
  )
ggsave('Extended_Figures/Ext_Fig3_Sinaplots_partners_nonpartners.pdf', height = 2, width = 2)
dev.off()

# make a table for the supplement showing the top correlations from the sinaplot
corr_for_sina %>% 
  filter(is_strong, partner) %>%
  select(mutant, strain, pearson, is_core) %>%
  arrange(desc(pearson)) %>% 
  rename('GSP1 mutant' = mutant,
         'Partner strain name' = strain,
         'Pearson Correlation Coefficient' = pearson,
         'Residue in core' = is_core) %>% 
  write_csv('Supplementary_Data_Tables/Supp_Table6_partner_correlations.csv')


