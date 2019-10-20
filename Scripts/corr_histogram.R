# histogram showing relationship between pval and FDR, for the supplement

library(tidyverse)
source('ucsf_colors.R')
# load('Data/spitzemapko_correlations_and_bonferroni_fdr_all.RData')


histdata <-
  correlations %>% 
  select(query_uniq1, query_uniq2, pearson, two_sided_p_value) %>% 
  mutate('sig' = ifelse(two_sided_p_value > 0.05, FALSE, TRUE))

ggplot(histdata, aes(x=pearson, fill = sig)) +
  geom_histogram(alpha = 0.3, bins = 5000) +
  scale_fill_manual(name = '', values = c(ucsf_colors$pink1, ucsf_colors$blue1), labels = c('p < 0.05', 'p >= 0.05')) +
  xlim(c(-0.25, 0.25)) + ggtitle('p-values of correlations of genetic interaction profiles') +
  xlab('Pearson Correlation Coefficient') + ylab('Count') +
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

ggsave('Supplemental_Figures/Correlations_histogram.png', height = 3, width = 3, dpi = 300)
dev.off()
