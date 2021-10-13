library(tidyverse)
library(ggforce)
source('ucsf_colors.R')

raw_data <-
  read_csv('Data/Westerns/westerns_data.csv', col_types = cols()) %>% 
  mutate('round' = factor(round), 'date' = factor(date))

# process the data
processed_data <-
  raw_data %>% 
  filter(to_include) %>% # remove some lanes and gels based on manual inspection
  filter(! sample %in% c('Blank','Ladder','T34N')) %>% 
  filter((gsp1 > 500000)  & (gsp1 < 3500000)) %>%  # keep Gsp1 within calibration thresholds
  group_by(date, gel) %>% 
  mutate('mean_WT_total_protein' = mean(total_protein[sample == 'WT'])) %>% 
  filter(!is.na(mean_WT_total_protein)) %>% # filter out 5 gels with WT too low (results in an NA here)
  mutate('total_protein_adjustment' = total_protein / mean_WT_total_protein) %>% 
  mutate('adjusted_gsp1' = gsp1/total_protein_adjustment) %>% 
  mutate('relative_expression' = adjusted_gsp1/mean(adjusted_gsp1[sample == 'WT'])) %>% 
  ungroup() %>%
  group_by(sample, round) %>%
  mutate('avg_rel_expr_by_round' = mean(relative_expression), 'n_in_round' = n()) %>%
  ungroup() %>% 
  group_by(sample) %>% 
  mutate('avg_rel_expr' = mean(avg_rel_expr_by_round),
         'std_rel_expr' = sd(avg_rel_expr_by_round),
         'n_sample' = n()) %>% 
  ungroup() %>% 
  mutate(residue = case_when(!sample %in% c('WT', 'MAT-a') ~ substr(sample, 1, nchar(sample)-1), T ~ sample),
         resn = case_when(sample == 'WT' ~ 0, sample == 'MAT-a' ~ 1,
                          !sample %in% c('WT', 'MAT-a') ~ as.numeric(substr(sample, 2, nchar(sample)-1)))) %>% 
  mutate(residue = factor(residue, levels = arrange(., resn) %>% pull(residue) %>% unique())) %>% 
  group_by(sample) %>% 
  mutate(n_rounds = n_distinct(round)) %>% 
  ungroup() %>% 
  mutate(std_rel_expr = ifelse(n_rounds > 2, std_rel_expr, NA))

# split mutants by strength for plotting purposes
strong_mutants <- c('H141E','Y157A','D79A','D79S','T34Q', 'T34E',
                    'K101R','T34G','T34A','R108L','R108I','Q147E',
                    'H141I','G80A','H141R','R112A','R112S','R108Y',
                    'R108Q','Y148I','R108G','R108A','R78K')
weak_mutants <- c('K143H', 'K169I', 'N105L', 'T139A', 'R108S', 'R108D',
                  'E115I', 'N105V', 'F58A', 'K129F', 'F58L', 'K129E',
                  'K129I', 'K154M', 'T34S', 'T34D', 'T34Y', 'A180T', 
                  'N84Y', 'T137G', 'Q147L', 'T139R', 'K129T', 'T34L',
                  'N102I', 'N102K', 'N102M', 'E115A', 'H141V', 'K143W',
                  'K143Y', 'K132H')
controls <- c('WT', 'MAT-a')


# for the figure, plot the strong mutants and weak mutants separately

strong_mutant_data <- filter(processed_data, sample %in% c(strong_mutants, controls))

ggplot(data = strong_mutant_data,
       aes(x = factor(sample),
           y = avg_rel_expr,
           fill = factor(ifelse(sample == 'WT', 'WT', 'mutant')))) +
  geom_col(position = 'dodge', color = 'black', width = 0.7) +
  geom_point(data = strong_mutant_data,
             aes(x = factor(sample),
                 y = avg_rel_expr_by_round),
             size = 0.5) +
  geom_errorbar(aes(ymin=avg_rel_expr-std_rel_expr,
                    ymax=avg_rel_expr+std_rel_expr),
                width = 0.4, position='dodge') +
  scale_fill_manual(values = c(ucsf_colors$gray3, ucsf_colors$navy1)) +
  geom_hline(yintercept=1, color = 'red', linetype = 'dashed') +
  xlab('S. cerevisiae strain with genomically integrated Gsp1 point mutant') + 
  ylab('Average Relative Expression MUT/WT') +
  ggtitle('Strong mutants by genetic interaction profiles') +
  facet_grid(~residue, scales = 'free_x', space = 'free_x', switch = 'x') +
  theme_classic() +
  theme(text=element_text(size=6, family='Helvetica'),
        axis.text.x=element_text(angle=45,hjust=1),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = 'none')
  
ggsave('Extended_Figures/Ext_Fig2A_Westerns_Strong_Mutants.pdf', width = 5, height = 2)
dev.off()

### print for Per Figure Source files
strong_mutant_data %>% write_tsv('Per_Figure_source_files/EDF2a.txt')
  
  
weak_mutant_data <- filter(processed_data, sample %in% c(weak_mutants, controls))

### print for Per Figure Source files
weak_mutant_data %>% write_tsv('Per_Figure_source_files/EDF2b.txt')


ggplot(data = weak_mutant_data,
       aes(x = factor(sample),
           y = avg_rel_expr,
           fill = factor(ifelse(sample == 'WT', 'WT', 'mutant')))) +
  geom_col(position = 'dodge', color = 'black', width = 0.7) +
  geom_point(data = weak_mutant_data,
             aes(x = factor(sample),
                 y = avg_rel_expr_by_round),
             size = 0.5) +
  geom_errorbar(aes(ymin=avg_rel_expr-std_rel_expr,
                    ymax=avg_rel_expr+std_rel_expr),
                width = 0.4, position='dodge') +
  scale_fill_manual(values = c(ucsf_colors$gray3, ucsf_colors$navy1)) +
  geom_hline(yintercept=1, color = 'red', linetype = 'dashed') +
  xlab('S. cerevisiae strain with genomically integrated Gsp1 point mutant') + 
  ylab('Average Relative Expression MUT/WT') +
  ggtitle('Weak mutants by genetic interaction profiles') +
  facet_grid(~residue, scales = 'free_x', space = 'free_x', switch = 'x') +
  theme_classic() +
  theme(text=element_text(size=6, family='Helvetica'),
        axis.text.x=element_text(angle=45,hjust=1),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = 'none')

ggsave('Extended_Figures/Ext_Fig2B_Westerns_Weak_Mutants.pdf', width = 5, height = 2)
dev.off()


# # plot all mutants
# processed_data %>%
#   ggplot(aes(x = factor(sample), y = avg_rel_expr,
#              fill = factor(ifelse(sample == 'WT', 'WT', 'mutant')))) +
#   geom_col(position = 'dodge', color = 'black', width = 0.8) +
#   geom_errorbar(aes(ymin=avg_rel_expr-std_rel_expr,
#                     ymax=avg_rel_expr+std_rel_expr),
#                 width = 0.4, position='dodge') +
#   scale_fill_manual(values = c(ucsf_colors$gray3, ucsf_colors$navy1)) +
#   geom_hline(yintercept=1, color = 'red', linetype = 'dashed') +
#   xlab('Gsp1 strain') + ylab('Average Relative Expression MUT/WT') +
#   facet_grid(~residue, scales = 'free_x', space = 'free_x', switch = 'x') +
#   theme_classic() +
#   theme(text=element_text(size=6, family='Helvetica'),
#         axis.text.x=element_text(angle=45,hjust=1),
#         strip.background = element_blank(),
#         strip.text.x = element_blank(),
#         legend.position = 'none')
# 
# ggsave('Supplemental_Figures/Westerns_All_Mutants.pdf', width = 14, height = 2)
# dev.off()


# Also plot westerns by residue, in one table, and one residue per page
# processed_data %>% 
#   mutate(id = paste(date, gel, lane, sep = '_')) %>% 
#   ggplot(aes(x = id, y = relative_expression, fill = round)) +
#   geom_col(position= 'dodge') +
#   facet_wrap(~sample, scales = 'free') +
#   ylab('Relative Expression MUT/WT') +
#   theme_classic() +
#   theme(text=element_text(size=6, family='Helvetica'),
#         axis.text.x=element_text(angle=90))
#ggsave('Supplemental_Figures/Westerns_by_residue.pdf', width = 8.5, height = 11)
#dev.off()

# plots <- list()
# for (i in seq(length(unique(processed_data$sample)))) {
#   sample_name <- unique(processed_data$sample)[i]
#   plots[[i]] <-
#     processed_data %>% 
#     mutate(id = paste(date, gel, lane, sep = '_')) %>% 
#     filter(sample == sample_name) %>% 
#     ggplot(aes(x = id, y = relative_expression, fill = round)) +
#     geom_col(position= 'dodge') +
#     xlab(sample_name) + ylab('Relative Expression MUT/WT') +
#     theme_classic() +
#     theme(axis.text.x=element_text(angle=90,hjust=1))
# }
# #pdf('Supplemental_Figures/Westerns_one_residue_per_page.pdf', onefile = TRUE)
# #plots
# #dev.off()


##### REVISIONS
# plot sina plots comparing weak to strong mutants across all mutants
# simple plot just comparing distribution
# processed_data %>%
#   mutate(mutant_strength = case_when(sample %in% strong_mutants ~ 'Strong Mutants',
#                                      sample %in% weak_mutants ~ 'Weak Mutants')) %>%
#   filter(!is.na(mutant_strength)) %>%
#   ggplot(aes(x = mutant_strength, y = relative_expression)) +
#   geom_sina() +
#   stat_summary(fun.y=mean, geom='point', color=ucsf_colors$pink1, size = 5)
# 
# ggsave('Revisions/Sina_westerns_all_points.pdf', height = 4, width = 4)
# dev.off()

# can also make a plot showing row means
data_for_sina <-
  processed_data %>% 
  mutate(mutant_strength = case_when(
    sample %in% strong_mutants ~ 'Strong Mutants',
    sample %in% weak_mutants ~ 'Weak Mutants')) %>% 
  filter(!is.na(mutant_strength))

# can also replace avg_rel_expr_by_round with relative_expression to plot every technical replicate
# (i.e. measurement from a single lane) as points
data_for_mean <-
  data_for_sina %>%
  select(mutant_strength, round, avg_rel_expr_by_round) %>%
  unique() %>%
  group_by(mutant_strength, round) %>% 
  mutate(avg_per_round_by_mutant_strength = mean(avg_rel_expr_by_round)) %>%
  select(-avg_rel_expr_by_round) %>%
  unique()
## get n
data_for_sina %>% group_by(mutant_strength) %>% 
  summarise('count' = n())

ggplot(data_for_sina,
       aes(x = mutant_strength,
       # y = relative_expression,
       y = avg_rel_expr_by_round,
       # color = round,
       group = mutant_strength)) +
  geom_sina(maxwidth=0.7, size=1.5, alpha = 0.75, stroke=0, seed = 2) +
  # scale_color_manual(name = 'Biological\nReplicate',
  #                    values = c(ucsf_colors$pink1, ucsf_colors$cyan1, ucsf_colors$green1)) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar",
               fatten = 0, width = 0.5, color = ucsf_colors$pink1, show.legend = FALSE) + 
  # geom_point(data = avg_per_round_by_mutant_strength,
  #            mapping = aes(mutant_strength, avg_per_round_by_mutant_strength,
  #                          fill = round, group=mutant_strength),
  #            pch=21, size=3,
  #            position = position_jitter(width = 0.02, seed = 3)) +
  # scale_fill_manual(name = 'Biological\nReplicate',
  #                   values = c(ucsf_colors$pink1, ucsf_colors$cyan1, ucsf_colors$green1),
  #                   ) +
  xlab('') + ylab('Relative Expression MUT/WT') +
  # guides(color = guide_legend(override.aes = list(size = 1.5))) +
  theme_classic() +
  theme(text=element_text(size=6, family='Helvetica'),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = 'top')

ggsave('Revisions/Extended_Figures/EDF_2/Ext_Fig2C_Sina_westerns_avg_per_bio_replicate_2.0.pdf',
       height = 4.1, width = 2, useDingbats=F)
dev.off()

### print for Per Figure Source file
data_for_sina %>% write_tsv('Per_Figure_source_files/EDF2c.txt')
