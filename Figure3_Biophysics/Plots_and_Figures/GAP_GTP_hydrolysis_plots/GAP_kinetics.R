library(tidyverse)
library(ggrepel)
directory <- 'GAP_GTP_hydrolysis_plots'
source('~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/ucsf_colors.R')

barplot_colors <- c(ucsf_colors$orange1, ucsf_colors$gray2, ucsf_colors$navy2)

MM.data <- read_tsv(file.path(directory, "GAP_kinetics_MichaelisMenten_parameters.txt"))
mut_ordered_by_kcat_Km <- MM.data %>% select(mutant, kcat_Km) %>% arrange(kcat_Km) %>% unique() %>% pull(mutant)
GAP_interface_mutations <- c('K132H')
YRB1_interface_mutations <- c('T34E', 'T34A', 'T34G', 'T34Q', 'T34L', 'T34S', 'A180T', 'F58A')
MM.data <- MM.data %>% 
  mutate("interface" = ifelse( (mutant %in% GAP_interface_mutations), "in RanGAP interface", "not in RanGAP interface" )) %>%
  mutate('interface' = ifelse( mutant == 'WT', 'WT', interface)) %>% 
  unique() %>% 
  mutate("mutant" = factor(mutant, mut_ordered_by_kcat_Km)) %>% 
  arrange(mutant)
MM.data

MM.data %>% 
  ggplot(aes(mutant, kcat_Km, fill = interface)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = kcat_Km - kcat_Km_sd, ymax = kcat_Km + kcat_Km_sd), width = 0.5, size = 0.2, alpha = 0.75) +
  ylab(expression("kcat / Km [s "^-{}^{1}*mu*"M"^-{}^{1}*"]")) +
  xlab('\npoint mutation in Ran') +
  ggtitle('RanGAP mediated Ran:GTP to Ran:GDP + Pi hydrolysis') + 
  labs(fill = element_blank()) +
  theme_classic() +
  scale_fill_manual(values = barplot_colors) +
  geom_hline(yintercept = MM.data$kcat_Km[MM.data$mutant == 'WT'], linetype = "dashed", size = 0.3) +
  theme(text = element_text(size = 16), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 18), 
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        axis.line = element_line(size = 0.1),
        legend.position = c(0.18, 0.82),
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "gray"))

ggsave(filename = file.path(directory, "GAP_kcat_over_Km.pdf"), height = 5, width = 9)


MM.data %>% 
  ggplot(aes(mutant, kcat_Km, fill = interface)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = kcat_Km - kcat_Km_sd, ymax = kcat_Km + kcat_Km_sd), width = 0.5, size = 0.2, alpha = 0.75) +
  coord_flip(clip = "off") +
  ylab(expression("kcat / Km [s "^-{}^{1}*mu*"M"^-{}^{1}*"]")) +
  xlab('\npoint mutation in Ran') +
  ggtitle('RanGAP mediated Ran:GTP hydrolysis') + 
  labs(fill = element_blank()) +
  theme_classic() +
  scale_fill_manual(values = barplot_colors) +
  geom_hline(yintercept = MM.data$kcat_Km[MM.data$mutant == 'WT'], linetype = "dashed", size = 0.3) +
  theme(text = element_text(size = 12), 
        legend.text = element_text(size = 12),
        axis.line = element_line(size = 0.1),
        legend.position = c(0.7, 0.2),
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "gray"))

ggsave(filename = file.path(directory, "horizontal_GAP_kcat_over_Km.pdf"), height = 8, width = 5)


MM.data %>% 
  ggplot(aes(mutant, kcat, fill = interface)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = kcat - kcat_sd, ymax = kcat + kcat_sd), width = 0.9) +
  ylab(expression("kcat [s "^-{}^{1}*"]")) +
  xlab('\npoint mutation in Ran') +
  ggtitle('RanGAP mediated Ran:GTP to Ran:GDP + Pi hydrolysis') + 
  labs(fill = element_blank()) +
  theme_classic() +
  scale_fill_manual(values = barplot_colors) +
  geom_hline(yintercept = MM.data$kcat[MM.data$mutant == 'WT'], linetype = "dashed", size = 0.3) +
  theme(text = element_text(size = 16), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 18), 
        legend.text = element_text(size = 15),
        axis.line = element_line(size = 0.1),
        legend.position = c(0.8, 0.85),
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "gray"))
ggsave(filename = file.path(directory, "GAP_kcat.pdf"), height = 5, width = 11)

MM.data %>% 
  ggplot(aes(mutant, Km, fill = interface)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = Km - Km_sd, ymax = Km + Km_sd), width = 0.5) +
  ylab(expression("Km ["*mu*"M]")) +
  xlab('\npoint mutation in Ran') +
  ggtitle('RanGAP mediated Ran:GTP to Ran:GDP + Pi hydrolysis') + 
  labs(fill = element_blank()) +
  theme_classic() +
  scale_fill_manual(values = barplot_colors) +
  geom_hline(yintercept = MM.data$Km[MM.data$mutant == 'WT'], linetype = "dashed", size = 0.3) +
  theme(text = element_text(size = 16), 
        axis.text.x = element_text(angle = 45, hjust = 1,  size = 18), 
        legend.text = element_text(size = 15),
        axis.line = element_line(size = 0.1),
        legend.position = c(0.2, 0.85),
        legend.background = element_rect(size = 0.1, linetype="solid", colour = "gray"))
ggsave(filename =  file.path(directory, "GAP_Km.pdf"), height = 5, width = 11)

MM.data %>% 
  ggplot(aes(x = kcat, y = log(Km),  color = interface)) + 
  geom_point(aes(size = kcat_Km_sd)) + 
  scale_color_manual(values = barplot_colors) +
  geom_text_repel(aes(label = mutant), show.legend = F, point.padding = 0.25) +
  theme_classic() +
  labs(color = '') +
  xlab(expression("kcat [s "^-{}^{1}*"]")) +
  ylab(label = expression("ln( Km ["*mu*"M] )")) +
  theme(text = element_text(size = 16),
        axis.line = element_line(size = 0.1)) +
  labs(size = 'std.dev(kcat/Km)') 
ggsave(filename = file.path(directory, "GAP_Km_kcat_scatterplot.pdf"), height = 5, width = 8)


