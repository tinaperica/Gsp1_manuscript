library(tidyverse)
library(viridis)
library(cowplot)
### This code makes supplementary file figures that show the GEF Michaelis Menten plots
### mutants with complete data
mutants <- c("WT","T34A", "T34E", "T34G", "T34L", "T34Q", "F58A", "R78K", 'D79A', "D79S", "G80A", "K101R",
             "R108A", "R108G", "R108I", "R108L", "R108Q", "R108Y", 'R112S', "K132H",
             "H141R","K143W", "Q147E", "Y157A", "A180T")
### GEF data - initial rates vs [Gsp1]
GEF_data <- read_tsv('Data/RanGEF_assay/GEF_assay_initial_rate_data.txt') %>% 
  separate(., sample, into = c('uniq', 'mutant'), sep = '_') %>% 
  mutate('date' = as.character(date)) %>% 
  filter(mutant %in% mutants)
max_v0 <- ceiling(max(GEF_data$v0))
## Michaelis Menten parameters
MM_parameters <- read_tsv('Data/RanGEF_assay/GEF_kinetics_MichaelisMenten_parameters.txt') %>% 
  filter(mutant %in% mutants) %>% 
  select(mutant, kcat, kcat_se, Km, Km_se, kcat_Km, kcat_Km_se)
#### make GEF MM example plots
plots <- list()
for (i in seq_along(mutants)) {
  mut <- mutants[i]
  temp <- GEF_data %>% 
    filter(mutant == mut)
  mm_param <- MM_parameters %>% filter(mutant == mut)
  predicted_data <- tibble('mutant' = mut, 'conc' = seq(0,1.1*max(temp$conc), 0.1), 'v0' = NA) %>% 
                mutate('predicted_mean' = (conc * mm_param$kcat) / (mm_param$Km + conc))
  plots[[mut]] <- temp %>% 
    mutate('predicted_mean' = (conc * mm_param$kcat) / (mm_param$Km + conc)) %>% 
    bind_rows(., predicted_data) %>% 
    mutate('predicted_plus' =  (conc * (mm_param$kcat + mm_param$kcat_se)) / ((mm_param$Km + mm_param$Km_se) + conc),
           'predicted_minus' =  (conc * (mm_param$kcat - mm_param$kcat_se)) / ((mm_param$Km - mm_param$Km_se) + conc)) %>% 
    ggplot(aes(x = conc, y = v0)) +
    geom_point(size = 2) +
    geom_line(aes(x = conc, y = predicted_mean), color = 'black') +
    geom_line(aes(x = conc, y = predicted_plus), color = 'gray') +
    geom_line(aes(x = conc, y = predicted_minus), color = 'gray') +
    theme_classic() +
    ylim(c(0, max_v0)) +
    xlab(expression('Gsp1 concentration  ['*mu*'M]')) +
    ylab(expression('v'[0]*" / GEF concentration in "*mu*"M  [s"^-{1}*"]")) +
    annotate("text", x = 1.2, y = 6, size = 2,
             label = "paste(k[cat], \" = \")", parse = T) +
    annotate("text", x = 4.3, y = 6, size = 2,
             label = str_c(round(mm_param$kcat, 1), " +/- ", round(mm_param$kcat_se, 2))) +
    annotate("text", x = 1.2, y = 5.6, size = 2,
             label = "paste(K[m], \" = \")", parse = T) +
    annotate("text", x = 4.3, y = 5.6, size = 2,
             label = str_c(round(mm_param$Km, 1), " +/- ", round(mm_param$Km_se, 2))) +
    theme(text = element_text(size = 7, family = 'Helvetica'), 
          axis.text.x = element_text(size = 7), 
          axis.title = element_text(size = 7),
          axis.line = element_line(size = 0.1),
          axis.ticks = element_line(size = 0.1)) +
    ggtitle(str_c(mut, ' Gsp1    n = ', count(temp)$n))
}
assembled_figure <- plot_grid(plotlist = plots, align = 'hv', nrow = 9, labels = 'auto', label_size = 8)
pdf('Revisions/Supplementary_Files/All_GEF_MichaelisMenten_plots.pdf',  width = 7, height = 24.75)
assembled_figure
dev.off()



