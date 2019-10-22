library(tidyverse)
library(viridis)
library(cowplot)
### plot WT, T34E and R108G data for GAP and GEF assay
### make supplementary figures

### GEF data
GEF_data <- read_tsv('Data/RanGEF_assay/fit_data.txt', guess_max = 60000) %>% 
  filter(sample %in% c('PE1_WT', 'PE3_T34E', 'PE26_R108Q')) %>% 
  select(date, sample, conc, GEF_conc, fit, Time, fluorescence, observed, predicted, condition, v0) %>% 
  separate(., sample, into = c('uniq', 'mutant'), sep = '_') %>% 
  mutate('group' = str_c(mutant, GEF_conc, fit, sep = ' ')) %>% 
  mutate('date' = as.character(date))
### first make raw data curves
choose_dates <- GEF_data %>% 
  select(date, mutant, conc, GEF_conc) %>% 
  unique() %>% 
  group_by(date, mutant) %>% 
  summarize('count' = n()) %>% 
  arrange(desc(count))
GEF_data_select <- GEF_data %>% 
  filter(date %in% c('20180827', '20181106', '20180913'))
groups <- GEF_data_select %>% pull(group) %>% unique()
mutants <- GEF_data_select %>% pull(mutant) %>% unique()
### make raw data example plots
plots <- list()
for (i in seq_along(mutants)) {
  mut <- mutants[i]
  to_plot <- GEF_data_select %>% filter(mutant == mut) 
  GEF_conc <- to_plot %>% pull(GEF_conc) %>% unique()
  rel_fluorescence <- to_plot %>% 
    group_by(condition) %>% 
    summarize('max' = max(fluorescence), 'min' = min(fluorescence), 'max_pred' = max(predicted), 'min_pred' = min(predicted))
  total_max <- max(rel_fluorescence$max)
  to_plot <- to_plot %>% inner_join(., rel_fluorescence, by = 'condition') %>% 
    mutate('relative_fluorescence' = (fluorescence-min)/(max-min), 'relative_predicted' = (predicted-min_pred)/(max_pred-min_pred)) %>% 
    select(condition, conc, fluorescence, relative_fluorescence, predicted, relative_predicted, Time, fit)
  plots[[mut]][['raw']] <- to_plot %>% 
    filter(Time < 2000) %>% 
    ggplot(aes(x = Time, y = fluorescence, color = conc)) + 
    geom_point(size = 0.08) +
    theme_classic() +
    xlab('Time / s') +
    ylab('Fluorescence') +
    labs(color = expression("Gsp1 / "*mu*"M")) +
    scale_color_viridis() +
    ggtitle(str_c('Gsp1 ', mut, '  ', GEF_conc, ' nM GEF')) +
    theme(text = element_text(size = 7, family = 'Helvetica'), 
                axis.text.x = element_text(size = 7), 
                axis.title = element_text(size = 7),
                legend.text = element_text(size = 6),
                legend.key.size =  unit(0.5, 'cm'),
                axis.line = element_line(size = 0.1),
                legend.position = 'bottom',
                legend.background = element_rect(size = 0.1, linetype = "solid", colour = "gray"))
  ###### plot some fits
  if (mut == 'R108L') {
    to_plot_fit <- to_plot %>% 
      filter(conc > 3 & conc < 10)
  } else {
    to_plot_fit <- to_plot %>% 
      filter(fit == 'exp' & conc %in% c(1, 2, 3))
  }
  plots[[mut]][['fit']] <- to_plot_fit %>% 
    filter(Time < 2000) %>% 
    ggplot(aes(x = Time, y = predicted, group = condition)) + 
    geom_line(size = 0.2) +
    geom_point(aes(x = Time, y = fluorescence, color = as.character(conc)), size = 0.3, alpha = 0.7) +
    theme_classic() +
    xlab('Time / s') +
    ylab('Fluorescence') +
    labs(color = expression("Gsp1 / "*mu*"M")) +
    scale_color_viridis(discrete = T) +
    ggtitle(str_c('Gsp1 ', mut, '   ', GEF_conc, ' nM GEF - fit data')) +
    theme(text = element_text(size = 7, family = 'Helvetica'), 
          axis.text.x = element_text(size = 7), 
          axis.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          legend.key.size =  unit(0.5, 'cm'),
          axis.line = element_line(size = 0.1),
          legend.position = 'bottom',
          legend.background = element_rect(size = 0.1, linetype = "solid", colour = "gray"))
}
###
assembled_figure <- plot_grid(plots$WT$raw, plots$WT$fit, 
                              plots$T34E$raw, plots$T34E$fit, 
                              plots$R108G$raw, plots$R108G$fit,
                              align = 'hv', nrow = 3, ncol = 2)

pdf('Supplemental_Figures/Supplemental_Fig8_example_of_raw_and_fit_GEF_data.pdf',  width = 8, height = 9)
assembled_figure
dev.off()

MM_parameters <- read_tsv('Data/RanGEF_assay/MM_data.txt') %>% 
  select(mutant, kcat, Km, kcat_sd, Km_sd) %>% 
  filter(mutant %in% c('WT', 'T34E', 'R108Q'))
GEF_data_MM <- GEF_data %>% 
  select(condition, mutant, conc, GEF_conc, v0) %>% 
  unique() %>% 
  select(mutant, conc, v0) %>% 
  arrange(mutant, conc)
#### make GEF MM example plots
plots <- list()
for (i in seq_along(mutants)) {
  mut <- mutants[i]
  temp <- GEF_data_MM %>% 
    filter(mutant == mut)
  mm_param <- MM_parameters %>% filter(mutant == mut)
  predicted_data <- tibble('mutant' = mut, 'conc' = seq(0,1.1*max(temp$conc), 0.1), 'v0' = NA) %>% 
                mutate('predicted_mean' = (conc * mm_param$kcat) / (mm_param$Km + conc))
  plots[[mut]] <- temp %>% 
    mutate('predicted_mean' = (conc * mm_param$kcat) / (mm_param$Km + conc)) %>% 
    bind_rows(., predicted_data) %>% 
    mutate('predicted_plus' =  (conc * (mm_param$kcat + mm_param$kcat_sd)) / ((mm_param$Km + mm_param$Km_sd) + conc),
           'predicted_minus' =  (conc * (mm_param$kcat - mm_param$kcat_sd)) / ((mm_param$Km - mm_param$Km_sd) + conc)) %>% 
    ggplot(aes(x = conc, y = v0)) +
    geom_point(size = 2) +
    geom_line(aes(x = conc, y = predicted_mean), color = 'black') +
    geom_line(aes(x = conc, y = predicted_plus), color = 'gray') +
    geom_line(aes(x = conc, y = predicted_minus), color = 'gray') +
    theme_classic() +
    ylim(c(0, 4)) +
    xlab(expression('Gsp1 concentration  ['*mu*'M]')) +
    ylab(expression('v'[0]*" / GEF concentration in "*mu*"M  [s"^-{1}*"]")) +
    theme(text = element_text(size = 7, family = 'Helvetica'), 
          axis.text.x = element_text(size = 7), 
          axis.title = element_text(size = 7),
          axis.line = element_line(size = 0.1),
          axis.ticks = element_line(size = 0.1)) +
    ggtitle(str_c(mut, ' Gsp1    n = ', count(temp)$n))
}
assembled_figure <- plot_grid(plots$WT, plots$T34E, plots$R108Q, align = 'h', nrow = 1, labels = 'auto', label_size = 8)
pdf('Supplemental_Figures/Supplemental_Fig9_GEF_MichaelisMenten_examples.pdf',  width = 7, height = 3)
assembled_figure
dev.off()



### GAP data
GAP_data <- read_tsv('Data/RanGAP_assay/all_GAP_data_parsed.txt') %>% 
  filter(sample %in% c('PE1_WT', 'PE3_T34E', 'PE4_R108L')) %>% 
  separate(., sample, into = c('uniq', 'mutant'), sep = '_') %>% 
  mutate('date' = as.character(date))
conditions <- GAP_data %>% select(mutant, condition) %>% unique() %>% arrange(desc(mutant))
discarded_conditions <- c('20190520-PE1_WT-C5-8-0-20', '20190520-PE1_WT-C6-8-1-20', '20190715-PE1_WT-A1-8-0-70', '20190715-PE1_WT-A2-8-1-70')
GAP_data <- GAP_data %>% 
  filter(! condition %in% discarded_conditions) %>% 
  mutate('experimental condition' = ifelse(GAP_conc != 0, str_c('    ', conc, '              ', GAP_conc), 'no GAP')) %>% 
  arrange(mutant, conc)
plots <- list()
for (i in seq_along(mutants)) {
  mut <- mutants[i]
  to_plot <- GAP_data %>% 
    filter(mutant == mut) %>% 
    filter(Time < 3000)
  no_GAP <- to_plot %>% 
    filter(GAP_conc == 0)
  with_GAP <- to_plot %>% 
    filter(GAP_conc > 0)
  plots[[mut]] <- with_GAP %>% 
    ggplot(aes(x = Time, y = corrected_fluorescence, color = `experimental condition`)) +
    geom_point(size = 0.3) +
    geom_point(data = no_GAP, aes(x = Time, y = corrected_fluorescence, color = `experimental condition`), size = 0.3) +
    theme_classic() +
    labs(color = expression("    Gsp1 / "*mu*"M   GAP / nM")) +
    xlab('Time / s') +
    ylab('Background corrected fluorescence') +
    scale_color_viridis(discrete = T) +
    ggtitle(str_c('Gsp1 ', mut)) +
    theme(text = element_text(size = 7, family = 'Helvetica'), 
          axis.text.x = element_text(size = 7), 
          axis.title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.key.size =  unit(0.4, 'cm'),
          axis.line = element_line(size = 0.1),
          legend.position = 'bottom',
          legend.direction = 'vertical',
          legend.background = element_rect(size = 0.1, linetype = "solid", colour = "gray"))
    
}

assembled_figure <- plot_grid(plots$WT, plots$T34E, plots$R108L, align = 'h', nrow = 1, labels = 'auto', label_size = 8)
pdf('Supplemental_Figures/Supplemental_Fig10_GAP_raw_data_examples.pdf',  width = 7, height = 4)
assembled_figure
dev.off()




