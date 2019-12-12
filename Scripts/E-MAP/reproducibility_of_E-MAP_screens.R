# make the E-MAP reproducibility plots
library(tidyverse)
library(pracma)
# this is a folder with all the individual screens (triplicates within one screen are merged)
screens_dir <- 'Data/E-MAP/individual_Gsp1_E-MAP_screens'
screens <- file.path(screens_dir, dir(screens_dir))
avg_merged_data <- read_tsv('Data/E-MAP/gsp1_pEMAP_avg_merged_gene_names.txt') %>% 
  gather('library_gene', 'avg_score', -Gene)
combined_screens <- tibble('Gene' = character(), 'library_gene' = character(), 'score' = double(), 'screen_date' = character(), 'avg_score' = double())
for (i in seq_along(screens)) {
  screen_date <- substr(screens[i], start = 42, stop = 49)
  # in the unmerged screens, some of the library genes are there twice, as internal control and they are merged in average merged files
  # when reading in, they are given unique names e.g. CCW12_1 - deal with it after gathering
  screen <- read_tsv(screens[i]) %>% 
    gather('library_gene', 'score', -Gene) %>% 
    mutate('screen_date' = screen_date) %>% 
    mutate('library_gene' = ifelse( substr(library_gene, nchar(library_gene)-1, nchar(library_gene)) == '_1', substr(library_gene, 1, nchar(library_gene)-2), library_gene)) %>%  
    inner_join(., avg_merged_data, by = c('Gene', 'library_gene'))
  combined_screens <- bind_rows(combined_screens, screen)
}
replicates <- combined_screens %>% 
  group_by(Gene, library_gene) %>% 
  summarise('replicate_count' = n()) %>% 
  ungroup() %>% 
  group_by(Gene) %>% 
  summarise('replicate_count' = min(replicate_count)) %>% 
  arrange(replicate_count)

  
to_fit <- combined_screens %>% 
  select(score, avg_score) %>% 
  filter(complete.cases(.))
odr_linear_fit <- odregress(to_fit$avg_score, to_fit$score)
odr_slope <- signif(odr_linear_fit$coeff[1], 2)
lm_fit <- lm(score ~ avg_score, data = to_fit)
adj.R.squared <- signif(summary(lm_fit)$adj.r.squared, 2)
pearson <- signif(cor(to_fit$avg_score, to_fit$score), 2)
combined_screens %>% 
  ggplot(aes(x = avg_score, y = score)) + geom_point(alpha = 0.4) +
  theme_classic() +
  xlab('final S-score - average score from 6-10 replicates') +
  ylab('genetic itneraction score from a single screen') +
  geom_abline(slope = odr_linear_fit$coeff[1], intercept = odr_linear_fit$coeff[2]) +
  annotate("text", x = -15, y = 8, label = str_c('R^2 = ', adj.R.squared, '\n', 'ODR fit slope = ', odr_slope, '\n', 'Pearson corr = ', pearson)) +
  theme(text = element_text(size = 7, family = 'Helvetica'), 
        axis.text = element_text(size = 7), 
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.2))

ggsave('Supplemental_Figures/Supp_Fig2_E-MAP_reproducibility_avg_score_vs_single_score.png', width = 4, height = 4)
