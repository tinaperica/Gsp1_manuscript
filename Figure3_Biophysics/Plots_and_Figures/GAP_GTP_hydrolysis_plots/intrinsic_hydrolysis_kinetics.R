library(tidyverse)
library(ggrepel)
directory <- 'GAP_GTP_hydrolysis_plots'
ucsf_colors <- list('orange1' = '#F48024', 'orange2' = '#F7A665', 'orange3' = '#FBCCA7', 
                    'cyan1' = '#18A3AC', 'cyan2' = '#5DBFC5', 'cyan3' = '#A3DADE', 
                    'navy1' = '#052049', 'navy2' = '#506380', 'navy3' = '#9BA6B6', 
                    'blue1' = '#178CCB', 'blue2' = '#5DAFDB', 'blue3' = '#A2D1EA',
                    'green1' = '#90BD31', 'green2' = '#B1D16F', 'green3' = '#D3E4AD', 
                    'purple1' = '#716FB2', 'purple2' = '#9C9AC9', 'purple3' = '#C6C5E0',
                    'yellow1' = '#FFDD00', 'yellow2' = '#FFE74D', 'yellow3' = '#FF199',
                    'pink1' = '#ED1848', 'pink2' = '#F25D7F', 'pink3' = '#F7A3B6',
                    'gray1' = '#4D4D4D', 'gray2' = '#999999', 'gray3' = '#B4B9BF')
barplot_colors <- c(ucsf_colors$gray2, ucsf_colors$navy2)
MM.data <- read_tsv(file.path(directory, "GAP_kinetics_MichaelisMenten_parameters.txt"))
intrinsic_hydrolysis <- read_tsv(file.path(directory, 'intrinsic_hydrolysis.txt')) %>% 
  inner_join(., MM.data, by = 'mutant') %>% 
  mutate('group' = ifelse( mutant == 'WT', 'WT', 'mutant')) %>% 
  arrange(mutant)

mut_ordered_by_intrinsic <- intrinsic_hydrolysis %>% 
  select(mutant, mean_rel_rate) %>% 
  arrange(mean_rel_rate) %>% 
  unique() %>% pull(mutant)

intrinsic_hydrolysis %>% 
  mutate("mutant" = factor(mutant, mut_ordered_by_intrinsic)) %>% 
  ggplot(aes(mutant, mean_rel_rate, fill = group)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = mean_rel_rate - sd_rel_rate, ymax = mean_rel_rate + sd_rel_rate), width = 0.5) +
  ylab(expression("intrinsic rate / [Ran:GTP] [s "^-{}^{1}*mu*"M]")) +
  xlab('\npoint mutation in Ran') +
  ggtitle('Intrinsic Ran:GTP to Ran:GDP + Pi hydrolysis\n') + 
  theme_light() +
  scale_fill_manual(values = barplot_colors) +
  geom_hline(yintercept = intrinsic_hydrolysis$mean_rel_rate[intrinsic_hydrolysis$mutant == 'WT'], linetype = "dashed", size = 0.3) +
  theme(text = element_text(size = 16), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = 'none')
ggsave(filename = file.path(directory, "intrinsic_hydrolysis_barplot.pdf"), height = 5, width = 9)

