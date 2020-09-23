#### plot point plot with fold change
library(tidyverse)
library(ggforce)
source('ucsf_colors.R')


data <- read_tsv('Data/APMS_data.txt') %>% 
  filter(! Prey_gene_name == 'GSP1') %>% 
  filter(norm == 'eqM') %>% 
  mutate('Prey_gene_name' = str_c(substr(Prey_gene_name, 1, 1),
                                  tolower(substr(Prey_gene_name, 2, nchar(Prey_gene_name))))) %>% 
  mutate('interface_partner' = str_c(substr(interface_partner, 1, 1),
                                     tolower(substr(interface_partner, 2, nchar(interface_partner))))) %>% 
  mutate('Prey_gene_name' = case_when(is.na(Prey_gene_name) ~ PreyORF,
                                      !is.na(Prey_gene_name) ~ Prey_gene_name)) %>% 
  mutate('sample' = ifelse(tag == 'N', str_c('N-3xFL-', mutant), str_c(mutant, '-C-3xFL'))) %>% 
  select(sample, tag, mutant, residue, Prey_gene_name, log2FC, adj.pvalue, deltarASA, interface_partner) %>% 
  unique() %>% 
  complete(interface_partner, nesting(sample, tag, mutant, residue, Prey_gene_name, log2FC, adj.pvalue), fill = list(deltarASA = 0))

core_partners <- data %>% pull(interface_partner) %>% unique()
ordered_mutants <- data %>% select(residue, mutant) %>% unique() %>% arrange(residue) %>% pull(mutant)
preys <- data %>% pull(Prey_gene_name) %>% unique()

# WHATS MISSING: LOG2FC = NA FOR PARTNERS THAT ARE NOT SEEN WITH THAT MUTANT

data_interfaces <-
  data %>% 
  filter(Prey_gene_name %in% core_partners & interface_partner %in% preys & Prey_gene_name == interface_partner) %>% 
  mutate(is_core = ifelse(deltarASA >= 0.25, T, F)) %>% 
  mutate(is_core = factor(is_core, levels = c(T, F)))

data_interfaces %>% 
  ggplot(aes(x = is_core, y = log2FC, color = is_core, size = adj.pvalue)) +
  geom_sina(alpha = 1, maxwidth = 0.5, seed = 2) +
  geom_hline(yintercept=0, linetype='dashed', color = "red") +
  scale_x_discrete(breaks=c(T, F), labels = c('mutant in\npartner interface',
                                              'mutant not in\npartner interface')) +
  scale_size(name = 'P-value', range = c(0.5, 0.01), breaks = c(1.0, 0.1, 0.01)) +
  scale_color_manual(name = 'category', guide = 'none',
                     values = c(ucsf_colors$gray1, ucsf_colors$gray3)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult=1), color = ucsf_colors$pink1,
               geom = "errorbar", width = 0.1, size = 0.75) +
  stat_summary(fun.y = mean, geom = "point", color = ucsf_colors$pink1, size = 2) +
  xlab('') +
  ylab('log2 fold change of\npartner pulled-down abundance') +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica", size = 6),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.ticks = element_line(size = 0.05),
    axis.ticks.length = unit(0.05, 'cm'),
    legend.position = 'right',
    axis.line = element_line(size = 0.1)
  )

ggsave('talks/APMS_sinaplot.pdf', height = 1.9, width = 3.1)
dev.off()

distr_interface <- filter(data_interfaces, is_core == T) %>% pull(log2FC)
distr_not_interface <- filter(data_interfaces, is_core == F) %>% pull(log2FC)

t.test(distr_interface, distr_not_interface)
