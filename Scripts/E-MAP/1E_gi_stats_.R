library(tidyverse)
library(ggforce)
source('ucsf_colors.R')

# load the spitzemapko, which contains the Gsp1 EMAP and the subset of the SGA
# that consists of genetic interactions where at least one strain is a knockout
load('Data/spitzemapko.rda')

# find intersecting arrays between the Gsp1 EMAP and the SGA
gsp1_arrays <-
  spitzemapko %>%
  filter(interaction_network == 'gsp1_pEMAP') %>%
  select(array_ORF) %>%
  distinct() %>%
  unlist()

sga_arrays <-
  spitzemapko %>%
  filter(interaction_network != 'gsp1_pEMAP') %>%
  select(array_ORF) %>%
  distinct() %>%
  unlist()

# Strong mutants are chosen from the Gsp1 EMAP clustering,
# split into two clusters
strong_mutants <- c('H141E','Y157A','D79A','D79S','T34Q', 'T34E',
                    'K101R','T34G','T34A','R108L','R108I','Q147E',
                    'H141I','G80A','H141R','R112A','R112S','R108Y',
                    'R108Q','Y148I','R108G','R108A','R78K')


# get confidence interval values using the spline fit to Collins data

# To compare Gsp1 mutants and SGA "query genes", consider overlapping "array" genes
data <-
  spitzemapko %>%
  filter(! query_mutant %in% c('NTER3XFLAG WT', 'CTER3XFLAG WT', 'T34N')) %>% 
  filter(array_ORF %in% intersect(gsp1_arrays, sga_arrays)) %>%
  mutate(group = case_when(  # Define a new variable "group" for plotting purposes
    interaction_network == 'ExN_SGA' ~ 'Essential\ngene KDs',
    interaction_network == 'NxN_SGA' ~ 'Non-essential\ngene KOs',
    (interaction_network == 'gsp1_pEMAP' & query_mutant %in% strong_mutants) ~ 'Strong\nGsp1 mutants',
    (interaction_network == 'gsp1_pEMAP' & ! query_mutant %in% strong_mutants) ~ 'Weak\nGsp1 mutants'))

data_filtered <- 
  data %>% 
  group_by(query_allele_name) %>%
  # mutate('n_sig' = sum((score < -4) | (score > 4.9))) %>% # Thresholds for significant interactions, 95% confidence
  mutate('n_sig' = sum((score < -3) | (score > 3))) %>% # Thresholds for significant interactions, slightly less stringent
  select(query_allele_name, n_sig, group) %>%
  unique()
  
data_filtered %>% 
  filter(grepl('GSP1', query_allele_name)) %>% 
  unique() %>% 
  arrange(desc(n_sig)) %>% 
  write_csv('Figure1_E-MAP/Num_sig_interactions_by_mutant_abs_score_threshold_3.csv')

gsp1_groups <- c('Strong\nGsp1 mutants', 'Weak\nGsp1 mutants')

# sina plots with total number of strong interactions per query
set.seed(2)
ggplot(data_filtered, aes(x = group, y = n_sig)) +
  geom_violin(data = filter(data_filtered, ! group %in% gsp1_groups),
              fill = ucsf_colors$gray3) +
  geom_jitter(data = filter(data_filtered, group %in% gsp1_groups),
              height = 0, width = 0.1, size = 1) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar",
               fatten = 0, width = 0.5, color = ucsf_colors$pink1) + 
  xlab('') +
  ylab('Number of significant\ngenetic interactions') +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica", size = 6),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.ticks = element_line(size = 0.05),
    axis.ticks.length = unit(0.05, 'cm'),
    legend.position = 'none',
    axis.line = element_line(size = 0.1)
  )
ggsave('Figure1_E-MAP/Plots/1E_GI_Stats_Violin.pdf', width = 2.5, height = 1.9)
ggsave('Figure1_E-MAP/Plots/1E_GI_Stats_Violin_for_EDF.pdf', width = 2.2, height = 1.5)

ggsave('Figure1_E-MAP/Plots/1E_GI_Stats_Violin_FinalFormatting.pdf', width = 3, height = 1.3)

### For per figure source file
data_filtered %>% 
  mutate(group = str_replace(group, '\n',' ')) %>% 
  write_csv('Per_Figure_source_files/Fig1D.csv')

#dev.off()

# cdf plot for all interactions
data %>%
  ggplot(aes(x = score, color = group)) +
  geom_step(stat="ecdf") +
  scale_color_manual(values=c(ucsf_colors$navy1,
                              ucsf_colors$cyan1,
                              ucsf_colors$purple1,
                              ucsf_colors$green1)) +
  xlab('S-score') + ylab('Fraction of data') + xlim(-5, 5) + ylim(0,1) +
  ggtitle('Cumulative density of genetic interaction scores') +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica", size = 6),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.ticks = element_line(size = 0.05),
    axis.ticks.length = unit(0.05, 'cm'),
    legend.title = element_blank(),
    axis.line = element_line(size = 0.1)
    # plot.title = element_text(hjust = 0.5, size = 4)
  )
ggsave('Supplemental_Figures/GI_statistics_CDF.pdf', width = 3, height = 2)
dev.off()
