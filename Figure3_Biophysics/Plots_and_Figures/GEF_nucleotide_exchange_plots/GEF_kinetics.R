library(tidyverse)
library(ggrepel)
directory <- 'GEF_nucleotide_exchange_plots'
source('~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/ucsf_colors.R')

barplot_colors <- c(ucsf_colors$cyan1, ucsf_colors$gray2, ucsf_colors$navy2)

MM.data <- read_tsv(file.path(directory, "GEF_kinetics_MichaelisMenten_parameters.txt"))
mut_ordered_by_kcat_Km <- MM.data %>% select(mutant, kcat_Km) %>% arrange(kcat_Km) %>% unique() %>% pull(mutant)
GEF_interface_mutations <- c('K101R', 'R108L', 'R108I', 'R108Y', 'R018A', 'N105L', 'R108Q', 'R108G')
MM.data <- MM.data %>% 
  mutate("GEF interface" = ifelse( (mutant %in% GEF_interface_mutations), "in RanGEF interface", "not in RanGEF interface" )) %>%
  mutate('GEF interface' = ifelse( mutant == 'WT', 'WT', `GEF interface`)) %>% 
  unique() %>% 
  mutate("mutant" = factor(mutant, mut_ordered_by_kcat_Km)) %>% 
  arrange(mutant)
MM.data

MM.data %>% 
  ggplot(aes(mutant, kcat_Km, fill = `GEF interface`)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = kcat_Km - kcat_Km_sd, ymax = kcat_Km + kcat_Km_sd), width = 0.5, size = 0.2, alpha = 0.75) +
  ylab(expression("kcat / Km [s "^-{}^{1}*mu*"M"^-{}^{1}*"]")) +
  xlab('\npoint mutation in Ran') +
  ggtitle('RanGEF mediated nucleotide exchange (Ran:GDP to Ran:mant-GTP)') + 
  labs(fill = element_blank()) +
  theme_classic() +
  scale_fill_manual(values = barplot_colors) +
  geom_hline(yintercept = MM.data$kcat_Km[MM.data$mutant == 'WT'], linetype = "dashed", size = 0.3) +
  theme(text = element_text(size = 16), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 18), 
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        legend.position = c(0.18, 0.65),
        axis.line = element_line(size = 0.1),
        legend.background = element_rect(size = 0.1, linetype="solid", colour = "gray"))

ggsave(filename = file.path(directory, "GEF_kcat_over_Km.pdf"), height = 5, width = 9)


MM.data %>% 
  ggplot(aes(mutant, kcat_Km, fill = `GEF interface`)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = kcat_Km - kcat_Km_sd, ymax = kcat_Km + kcat_Km_sd), width = 0.5, size = 0.2, alpha = 0.75) +
  coord_flip(clip = "off") +
  ylab(expression("kcat / Km [s "^-{}^{1}*mu*"M"^-{}^{1}*"]")) +
  xlab('\npoint mutation in Ran') +
  ggtitle('RanGEF mediated nucleotide exchange (Ran:GDP to Ran:mant-GTP)') + 
  labs(fill = element_blank()) +
  theme_classic() +
  scale_fill_manual(values = barplot_colors) +
  geom_hline(yintercept = MM.data$kcat_Km[MM.data$mutant == 'WT'], linetype = "dashed", size = 0.3) +
  theme(text = element_text(size = 12), 
        legend.text = element_text(size = 12),
        axis.line = element_line(size = 0.1),
        legend.position = c(0.7, 0.2),
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "gray"))

ggsave(filename = file.path(directory, "horizontal_GEF_kcat_over_Km.pdf"), height = 8, width = 5)



MM.data %>% 
  ggplot(aes(mutant, kcat, fill = `GEF interface`)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = kcat - kcat_sd, ymax = kcat + kcat_sd), width = 0.5) +
  ylab(expression("kcat [s "^-{}^{1}*"]")) +
  xlab('\npoint mutation in Ran') +
  ggtitle('RanGEF mediated nucleotide exchange (Ran:GDP to Ran:mant-GTP)') + 
  labs(fill = element_blank()) +
  theme_light() +
  scale_fill_manual(values = barplot_colors) +
  geom_hline(yintercept = MM.data$kcat[MM.data$mutant == 'WT'], linetype = "dashed", size = 0.3) +
  theme(text = element_text(size = 16), 
        axis.text.x = element_text(angle = 45, hjust = 1,  size = 18), 
        legend.text = element_text(size = 9),
        legend.position = c(0.8, 0.85),
        axis.line = element_line(size = 0.1),
        legend.background = element_rect(size = 0.1, linetype = 'solid', colour = 'gray'))
ggsave(filename = file.path(directory, "GEF_kcat.pdf"), height = 5, width = 11)

all_Km <- MM.data %>% 
  ggplot(aes(mutant, Km, fill = `GEF interface`)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = Km - Km_sd, ymax = Km + Km_sd), width = 0.5) +
  ylab(expression("Km ["*mu*"M]")) +
  xlab('\npoint mutation in Ran') +
  ggtitle('RanGEF mediated nucleotide exchange (Ran:GDP to Ran:mant-GTP)') + 
  labs(fill = element_blank()) +
  theme_light() +
  scale_fill_manual(values = barplot_colors) +
  geom_hline(yintercept = MM.data$Km[MM.data$mutant == 'WT'], linetype = "dashed", size = 0.3) +
  theme(text = element_text(size = 16), 
        axis.text.x = element_text(angle = 45, hjust = 1,  size = 18), 
        legend.text = element_text(size = 9),
        legend.position = c(0.2, 0.85),
        axis.line = element_line(size = 0.1),
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "gray"))
ggsave(filename =  file.path(directory, "GEF_Km.pdf"), height = 5, width = 11)

cut_mut_ordered_by_kcat_Km <- mut_ordered_by_kcat_Km[4:27]
zoom_Km <- MM.data %>% 
  filter(mutant %in% cut_mut_ordered_by_kcat_Km) %>% 
  ggplot(aes(mutant, Km, fill = `GEF interface`)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = Km - Km_sd, ymax = Km + Km_sd), width = 0.5) +
  ylab(expression("Km ["*mu*"M]")) +
  xlab('\npoint mutation in Ran') +
  labs(fill = element_blank()) +
  theme_light() +
  scale_fill_manual(values = barplot_colors) +
  theme(text = element_text(size = 16), 
        axis.text.x = element_text(angle = 45, hjust = 1,  size = 18), 
        legend.text = element_text(size = 9),
        axis.line = element_line(size = 0.1),
        legend.position = "none",
        legend.background = element_rect(size = 0.1, linetype="solid", colour = "gray"))
ggsave(filename = file.path(directory, "GEF_Km_zoom.pdf"), width = 11, height = 5)





MM.data %>% 
  ggplot(aes(x = kcat, y = log(Km),  color = `GEF interface`)) + 
  geom_point(aes(size = kcat_Km_sd)) + 
  #geom_text(aes(label = mutant), show.legend = F) + 
  scale_color_manual(values = barplot_colors) +
  geom_text_repel(aes(label = mutant), show.legend = F, point.padding = 0.25) +
  theme_classic() +
  theme(text = element_text(size = 16),
        axis.line = element_line(size = 0.1)) +
  labs(color = '') +
  xlab(expression("kcat [s "^-{}^{1}*"]")) +
  ylab(label = expression("ln( Km ["*mu*"M] )")) +
  labs(size = 'std.dev(kcat/Km)') 
ggsave(filename = file.path(directory, "GEF_Km_kcat_scatterplot.pdf"), height = 5, width = 8)

