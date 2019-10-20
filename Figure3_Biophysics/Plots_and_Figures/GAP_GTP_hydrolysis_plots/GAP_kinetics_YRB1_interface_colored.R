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

barplot_colors <- c(ucsf_colors$gray2, ucsf_colors$yellow1, ucsf_colors$navy2)

MM.data <- read_tsv(file.path(directory, "GAP_kinetics_MichaelisMenten_parameters.txt"))
mut_ordered_by_kcat_Km <- MM.data %>% select(mutant, kcat_Km) %>% arrange(kcat_Km) %>% unique() %>% pull(mutant)
YRB1_interface_mutations <- c('T34E', 'T34A', 'T34G', 'T34Q', 'T34L', 'T34S', 'A180T', 'F58A')


MM.data <- MM.data %>% 
  mutate("interface" = ifelse( (mutant %in% YRB1_interface_mutations), "residue in interface\nwith YRB1", "not in YRB1 interface\n" )) %>%
  mutate('interface' = ifelse( mutant == 'WT', 'WT', interface)) %>% 
  unique() %>% 
  mutate("mutant" = factor(mutant, mut_ordered_by_kcat_Km)) %>% 
  arrange(mutant)
MM.data

MM.data %>% 
  ggplot(aes(mutant, kcat_Km, fill = interface)) + geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = kcat_Km - kcat_Km_sd, ymax = kcat_Km + kcat_Km_sd), width = 0.5) +
  ylab(expression("kcat / Km [s "^-{}^{1}*mu*"M]")) +
  xlab('\npoint mutation in Ran') +
  ggtitle('RanGAP mediated Ran:GTP to Ran:GDP + Pi hydrolysis\n') + 
  labs(fill = element_blank()) +
  theme_light() +
  scale_fill_manual(values = barplot_colors) +
  geom_hline(yintercept = MM.data$kcat_Km[MM.data$mutant == 'WT'], linetype = "dashed", size = 0.3) +
  theme(text = element_text(size = 16), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.text = element_text(size = 9),
        legend.position = c(0.12, 0.82),
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "gray"))

ggsave(filename = file.path(directory, "GAP_kcat_over_Km_YRB1_interface_mutants.pdf"), height = 5, width = 9)