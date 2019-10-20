library(tidyverse)
library(cowplot)

files <- dir(path = 'Data/C18_column_HPLC_data', pattern = '*.csv', full.names = T)
names(files) <- basename(files)

data <-
  files %>% 
  map_dfr(read.csv, .id = "filename") %>%
  as_tibble() %>% 
  mutate(sample = substr(filename, 1, nchar(filename)-4) %>%
           gsub(patter = '_', replacement = ' '))

plot_spectrum <- function(sample_name) {
  data %>% 
  filter(sample == sample_name, retention_minutes < 25) %>% 
  ggplot(aes(x = retention_minutes, y = abs_254nm)) +
    geom_line(size = 0.2) +
    ylab('254 nm absorbance') +
    ggtitle(sample_name) + 
    theme_classic() +
    theme(
      text = element_text(family = "Helvetica", size = 6),
      plot.title = element_text(hjust=1),
      axis.title = element_text(size = 5),
      axis.text = element_text(size = 5),
      axis.ticks = element_line(size = 0.05),
      axis.ticks.length = unit(0.05, 'cm'),
      axis.line = element_line(size = 0.1)
    )
}

p1 <-
  plot_spectrum('GTP GDP Mix 25uM') +
    xlab(NULL) +
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
p2 <-
  plot_spectrum('Gsp1 GTP WT 25uM') +
    xlab("Retention (minutes)")

plot_grid(p1, p2, nrow = 2, ncol = 1)

ggsave('Biophysics_Figure3/Biophysics_Figure_Panels/Biophysics_supplementary_figure_panels/HPLC_GTP_Loading.pdf',
       height = 2, width = 2)
