library(bio3d)
library(tidyverse)
library(RColorBrewer)

standard_SASA <- read_tsv("per_AA_standard_ASA_CWilke.txt", col_names = T)
switch_loops <- tibble(resi = c(19:25, 41:47, 69:77), "loop" = c(rep("ploop", 7), rep("switch1", 7), rep("switch2",9)))
pdb_dir <- "pdbs/clean"
chainA_dir <- "pdbs/clean_chainA/"
other_chain_dir <- "pdbs/clean_other_chain/"
files_list <- list.files(path = pdb_dir)
all_chains_file_paths <- file.path(pdb_dir, files_list)
chainA_file_paths <- file.path(chainA_dir, files_list)
other_chain_file_paths <- file.path(other_chain_dir, files_list)
Gsp1_sequence <- 'MSAPAANGEVPTFKLVLVGDGGTGKTTFVKRHLTGEFEKKYIATIGVEVHPLSFYTNFGEIKFDVWDTAGQEKFGGLRDGYYINAQCAIIMFDVTSRITYKNVPNWHRDLVRVCENIPIVLCGNKVDVKERKVKAKTITFHRKKNLQYYDISAKSNYNFEKPFLWLARKLAGNPQLEFVASPALAPPEVQVDEQLMQQYQQEMEQATALPLPDEDDADL'
index <- read_tsv("index.txt", col_names = T)
interface_SASA_table <- tibble(file_path = character(), chain = character(), resno = integer(), resid = character(),
              #SASA_difference = double(), percent_SASA_change = double(), 
              interface = character(), rosetta_num = integer())

for (i in seq_along(all_chains_file_paths)) {
  file_path <- all_chains_file_paths[i]
  chainA_path <- chainA_file_paths[i]
  other_chain_path <- other_chain_file_paths[i]
  pdb <- read.pdb(file_path, rm.alt = FALSE)
  pdbA <- read.pdb(chainA_path, rm.alt = FALSE)
  pdbB <- read.pdb(other_chain_path, rm.alt = FALSE)
  SASA_A_monomer <- dssp(pdbA)$acc
  SASA_B_monomer <- dssp(pdbB)$acc
  SASA_A_complex <- dssp(pdb)$acc[1:length(SASA_A_monomer)]
  SASA_B_complex <- dssp(pdb)$acc[(length(SASA_A_monomer) + 1):(length(SASA_A_monomer)+length(SASA_B_monomer))]
  coord.temp <- as_tibble(pdb$atom[, c("chain", "resno", "resid", "elety", "x", "y", "z")]) %>%
    filter(resid != "HOH" & resid != "SO4" & resid != "MG" & 
             resid != "GNP" & resid != "GDP"  & resid != "GTP" & resid != "AF3" & resid != "GOL" &
             resid != "EDO" & resid != "CL")
  chainA_deltaSASA <- coord.temp %>% 
    filter(chain == "A") %>% 
    select(chain, resno, resid) %>% unique() %>% 
    add_column(., SASA_A_monomer, SASA_A_complex, file_path) %>%
    inner_join(., standard_SASA, by = c("resid" = "AA")) %>% 
    select(-A, -Residue) %>% 
    mutate("rASAm" = SASA_A_monomer/maxASA, "rASAc" = SASA_A_complex/maxASA) %>% 
    mutate("deltarASA" = rASAm - rASAc) %>% 
    mutate("interface" = ifelse((rASAm >= 0.25 & rASAc < 0.25 & deltarASA >= 0.25), "core", "rim")) %>% 
    mutate("interface" = ifelse(rASAm < 0.25, "support", interface)) %>% 
    mutate('interface' = ifelse(deltarASA == 0, 'not_interface', interface)) %>% 
    select(file_path, chain, resno, resid, deltarASA, interface) %>% 
    bind_cols(., unique(tibble('aa' = pdbseq(pdbA), 'pdb_seq_number' = names(pdbseq(pdbA))))) %>% 
    filter(resno > 0)
  length <- length(chainA_deltaSASA$aa)
  chainA_deltaSASA <- chainA_deltaSASA %>% mutate('rosetta_num' = seq(1, length, 1))
  chainB_deltaSASA <- coord.temp %>% 
    filter(chain != "A") %>% 
    select(chain, resno, resid) %>% unique() %>% 
    add_column(., SASA_B_monomer, SASA_B_complex, file_path) %>%
    inner_join(., standard_SASA, by = c("resid" = "AA")) %>% 
    select(-A, -Residue) %>% 
    mutate("rASAm" = SASA_B_monomer/maxASA, "rASAc" = SASA_B_complex/maxASA) %>% 
    mutate("deltarASA" = rASAm - rASAc) %>% 
    mutate("interface" = ifelse((rASAm >= 0.25 & rASAc < 0.25 & deltarASA >= 0.25), "core", "rim")) %>% 
    mutate("interface" = ifelse(rASAm < 0.25, "support", interface)) %>% 
    mutate('interface' = ifelse(deltarASA == 0, 'not_interface', interface)) %>% 
    select(file_path, chain, resno, resid, deltarASA, interface) %>% 
    bind_cols(., unique(tibble('aa' = pdbseq(pdbB), 'pdb_seq_number' = names(pdbseq(pdbB))))) %>% 
    filter(resno > 0)
  length <- length(chainB_deltaSASA$aa)
  chainB_deltaSASA <- chainB_deltaSASA %>% mutate('rosetta_num' = seq(1, length, 1))
  interface_SASA_table <- interface_SASA_table %>% 
    bind_rows(., select(chainA_deltaSASA, -pdb_seq_number)) %>% 
    bind_rows(., select(chainB_deltaSASA, -pdb_seq_number))
}

interface_SASA_table <- index %>%
  select(file_path, chain, partner) %>% 
  full_join(., interface_SASA_table, by = c('file_path', 'chain')) %>% 
  mutate('partner' = ifelse(is.na(partner), 'GSP1', partner))
files <- interface_SASA_table %>% pull(file_path) %>% unique()

#contacts_table <- as_tibble(contacts_table)
#write_delim(contacts_table, "4A_contacts.txt", delim = "\t")

alignments_with_yeast_tibble <- tibble('file_path' = character(), 'chain' = character(), 'pdb_seq' = character(), 'yeast_seq' = character()) 
for (i in seq_along(files)) {
  file_path_temp <- files[i]
  seq1 <- interface_SASA_table %>% 
    filter(file_path == file_path_temp & chain != 'A') %>% pull(aa) %>% str_c(collapse = '')
  partner_chain <- interface_SASA_table %>% 
    filter(file_path == file_path_temp & chain != 'A') %>% pull(chain) %>% unique()
  index_temp <- index %>% 
    filter(file_path == file_path_temp)
  sequences <- seqbind(seq1, index_temp$`partner yeast sequence`)
  full_seq_aln <- seqaln(sequences) 
  aln_tibble <- tibble('file_path' = file_path_temp, 'chain' = partner_chain, 'pdb_seq' = full_seq_aln$ali[1,], 'yeast_seq' = full_seq_aln$ali[2,]) %>% 
    filter(! pdb_seq == '-')
  length <- length(aln_tibble$pdb_seq)
  aln_tibble <- aln_tibble %>% 
    mutate('rosetta_num' = seq(1, length, 1))
  alignments_with_yeast_tibble <- alignments_with_yeast_tibble %>% 
    bind_rows(., aln_tibble)
  seq1 <- interface_SASA_table %>% 
    filter(file_path == file_path_temp & chain == 'A') %>% pull(aa) %>% str_c(collapse = '')
  sequences <- seqbind(seq1, Gsp1_sequence)
  full_seq_aln <- seqaln(sequences) 
  aln_tibble <- tibble('file_path' = file_path_temp, 'chain' = 'A', 'pdb_seq' = full_seq_aln$ali[1,], 'yeast_seq' = full_seq_aln$ali[2,]) %>% 
    filter(! pdb_seq == '-')
  length <- length(aln_tibble$pdb_seq)
  aln_tibble <- aln_tibble %>% 
    mutate('rosetta_num' = seq(1, length, 1))
  alignments_with_yeast_tibble <- alignments_with_yeast_tibble %>% 
    bind_rows(., aln_tibble)
}
alignments_with_yeast_tibble <- index %>%
  select(file_path, chain, partner) %>% 
  full_join(., alignments_with_yeast_tibble, by = c('file_path', 'chain')) %>% 
  mutate('partner' = ifelse(is.na(partner), 'GSP1', partner))

merged <- interface_SASA_table %>% 
  inner_join(., alignments_with_yeast_tibble, by = c('file_path', 'chain', 'rosetta_num', 'partner')) %>% 
  mutate('identical' = ifelse(pdb_seq == yeast_seq, TRUE, FALSE))

interface_res <- merged %>% 
  filter(! interface == 'not_interface') %>% 
  group_by(file_path, chain, partner) %>% 
  summarize('count' = n())

interface_res_conserved <- merged %>% 
  filter(! interface == 'not_interface' & identical == T) %>% 
  group_by(file_path, chain, partner) %>% 
  summarize('conserved_n' = n())

merged_cons <- interface_res %>% 
  inner_join(., interface_res_conserved, by = c('file_path', 'chain', 'partner')) %>% 
  mutate('percent' = conserved_n/count*100)


merged %>% filter(grepl('3m1i1', file_path))
interface_SASA_table %>% filter(grepl('3m1i1', file_path))
