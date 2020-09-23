library(tidyverse)
interfaces <- read_tsv('Data/Gsp1_interfaces_SASA_and_conservation.txt')
index <- read_tsv('Scripts/complex_structure_analyses/index.txt')
interface_res <- interfaces %>% 
  filter(! interface == 'not_interface') %>% 
  group_by(file_path, chain, protein) %>% 
  summarize('count' = n())

interface_res_conserved <- interfaces %>% 
  filter(! interface == 'not_interface' & identical == T) %>% 
  group_by(file_path, chain, protein) %>% 
  summarize('conserved_n' = n())

merged_cons <- interface_res %>% 
  inner_join(., interface_res_conserved, by = c('file_path', 'chain', 'protein')) %>% 
  mutate('percent' = conserved_n/count*100)
write_tsv(merged_cons, path = 'Revisions/Supplementary_Files/Supplementary_File_1/interface_conservation.txt')

#### show that RNA1 interface is more conserved than the other interfaces
non_yeast_partners <- interfaces %>% 
  group_by(protein, chain, identical) %>% 
  summarise('count' = n()) %>% 
  filter(identical == FALSE & count > 1) %>% 
  pull(protein)
interfaces %>% 
  filter(protein != 'GSP1') %>% 
  filter(interface != 'not_interface' & protein %in% non_yeast_partners) %>% 
  mutate('residue' = ifelse(identical == T, 'conserved', 'not conserved')) %>% 
  ggplot(mapping = aes(y = deltaASA, x = rASAc, color = residue)) +
  geom_point(alpha = 0.5) +
  theme_classic() + 
  facet_wrap(~protein, nrow = 3) +
  ylab('surface area buried upon interface formation\n(ASAmonomer - ASAcomplex) / A2') +
  xlab('\nrelative accessible surface area in the complex') +
  theme(text = element_text(size = 7, family = 'Helvetica'), 
        axis.text.x = element_text(size = 7), 
        axis.text.y = element_text(size = 7), 
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size =  unit(0.4, 'cm'),
        axis.line = element_line(size = 0.1),
        axis.ticks = element_line(size = 0.1),
        legend.position = 'bottom',
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "gray"))
ggsave('Revisions/Supplementary_Files/Supplementary_File_1/Supplementary_File_1_Figures/conservation_of_interfaces_scatterplot.pdf', height = 5, width = 6)

#### count number of different interfaces in which each Gsp1 residue is part of the interface core
interface_core_count <- interfaces %>% 
  filter(file_path %in% index$file_path & protein == 'GSP1' & interface == 'core') %>% 
  group_by(yeast_num) %>% 
  summarise('number_of_interfaces' = n()) %>% 
  arrange(desc(number_of_interfaces))
Ran_Gsp1_index <- read_tsv('Scripts/complex_structure_analyses/Gsp1_to_Ran.index') %>% 
  filter(pdb == 'human_Ran') %>% select(-pdb) %>% rename(mamm_num = 'pdb_seq_num')
interface_core_count <- interface_core_count %>% inner_join(., Ran_Gsp1_index, by = c('yeast_num' = 'ref_seq_num'))
### categories of residues for pymol
### in zero interfaces
seq(1, 219, 1) %>% .[! . %in% interface_core_count$mamm_num] %>% str_c(., collapse = '+')
### in 1 interface
interface_core_count %>% filter(number_of_interfaces == 1) %>% pull(mamm_num) %>% str_c(., collapse = '+') 
### in 2-4 interfaces
interface_core_count %>% filter(number_of_interfaces > 1 & number_of_interfaces < 5) %>% pull(mamm_num) %>% str_c(., collapse = '+')
### in 5 or more interfaces
interface_core_count %>% filter(number_of_interfaces >= 5) %>% pull(mamm_num) %>% str_c(., collapse = '+')





interface_count <- interfaces %>% 
  filter(file_path %in% index$file_path & protein == 'GSP1' & interface != 'not_interface') %>% 
  group_by(yeast_num) %>% 
  summarise('number_of_interfaces' = n()) %>% 
  arrange(yeast_num)
interface_count <- interface_count %>% inner_join(., Ran_Gsp1_index, by = c('yeast_num' = 'ref_seq_num'))
### categories of residues for pymol
### in zero interfaces
seq(1, 219, 1) %>% .[! . %in% interface_count$mamm_num] %>% str_c(., collapse = '+')
### in 1 interface
interface_count %>% 
  filter(number_of_interfaces == 1) %>% 
  pull(mamm_num) %>% sort() %>% str_c(., collapse = '+') 
### in 2-4 interfaces
interface_count %>% 
  filter(number_of_interfaces > 1 & number_of_interfaces < 5) %>% 
  pull(mamm_num) %>% sort() %>% str_c(., collapse = '+')
### in 5 or more interfaces
interface_count %>% filter(number_of_interfaces >= 5) %>% 
  pull(mamm_num) %>% sort() %>% 
  str_c(., collapse = '+')

Gsp1_interfaces <- interfaces %>% 
  filter(protein == 'GSP1' & interface != 'not_interface') %>% 
  mutate('interface' = ifelse(interface != 'core', 'other', interface)) %>% 
  inner_join(., Ran_Gsp1_index, by = c('yeast_num' = 'ref_seq_num'))
#### show the residues in 3m1i
pymol_select <- function(part, partner_group, partner_count = 1, int = 'core') {
  selection_start_string <- str_c('select ', partner_group, '_', int, '_', partner_count, ', Gsp1 and resi') 
  temp <- Gsp1_interfaces %>% filter(partner %in% part & interface == int)
  line <- temp %>% 
    group_by(mamm_num) %>% 
    summarise('count' = n()) %>% 
    arrange(desc(count)) %>% 
    filter(count == partner_count) %>% 
    pull(mamm_num) %>% sort() %>% 
    str_c(., collapse = '+') %>% 
    append(selection_start_string, after = 0) %>% 
    str_c(., collapse = ' ')
  return(line)
}
## residues in Srm1
## residues in Rna1
## residues in Ntf2
## residues in Nup1/60
## residues in Yrb1 and Yrb2
## residues in Srp1
## residues in Kap95/Crm1/Los1/Kap104/Msn5/Cse1/Mtr10
residue_sets <- list(
  'GEF' = 'SRM1',
  'GAP' = 'RNA1',
  'Ntf2' = 'NTF2',
  'Nup1/60' = c('NUP1', 'NUP60'),
  'Yrb1/2' = c('YRB1', 'YRB2'),
  'Srp1' = 'SRP1',
  'karyopherins' = c('KAP95', 'CRM1', 'LOS1', 'KAP104', 'MSN5', 'CSE1', 'MTR10')
)
final_pymol_commands <- tibble('command' = character())
### if it's only one partner in the group, report core and other, if it's more than one partner in the group, 
## report only core residues, and group them by number of interfaces in which those residues are the core
for (partner_group in names(residue_sets)) {
  if (length(residue_sets[[partner_group]]) == 1) {
    partner <- residue_sets[[partner_group]]
    command <- pymol_select(part = partner, partner_group = partner_group, partner_count = 1, int = 'core')
    final_pymol_commands <- add_row(final_pymol_commands, command)
    command <- pymol_select(part = partner, partner_group = partner_group, partner_count = 1, int = 'other')
    final_pymol_commands <- add_row(final_pymol_commands, command)
  } else {
    n_partners <- length(residue_sets[[partner_group]])
    partners <- residue_sets[[partner_group]]
    for (i in seq(1, n_partners, 1)) {
      command <- pymol_select(part = partners, partner_group = partner_group, partner_count = i, int = 'core')
      final_pymol_commands <- add_row(final_pymol_commands, command)
    }
  }
}

write_tsv(final_pymol_commands, path = 'PyMOL_figures/Ran_interface_patches/interface_patches.txt')

mutations <- c(34, 58, 78, 79, 80, 84, 101, 102, 105, 108, 112, 115, 129, 132, 
              137, 139, 141, 143, 147, 148, 154, 157, 169, 180)
interface_cores <- interfaces %>% 
  filter(interface == 'core' & protein == 'GSP1' & yeast_num %in% mutations) %>% 
  select(partner, yeast_num, resno, aa, yeast_seq) %>% 
  arrange(yeast_num)

cons <- interfaces %>% 
  filter(protein == 'GSP1' & yeast_num %in% mutations) %>% 
  select(yeast_num, identical) %>% 
  arrange(yeast_num) %>% 
  unique()

  
### Make a table for Supplementary File 1 Table 2 
### (dSASA and core/rim/support call for each mutated residue)

table2 <- interfaces %>% 
  filter(protein == 'GSP1' & yeast_num %in% mutations) %>% 
  arrange(partner, yeast_num) %>% 
  select(partner, interface, deltarASA, yeast_num) %>% 
  mutate('int' = ifelse(deltarASA > 0,
              str_c(interface, round(deltarASA, 2), sep = ' / '),
              '')) %>% 
  select(partner, int, yeast_num) %>% 
  pivot_wider(names_from = partner, values_from = int)
write_tsv(table2, 'Revisions/Supplementary_Files/Supplementary_File_1/excel files/table2.txt')  

                        