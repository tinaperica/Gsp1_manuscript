library(bio3d)
library(tidyverse)

#### standard SAS values for each amino acid - needed to calculate relative ASA
standard_SASA <- read_tsv("Scripts/complex_structure_analyses/per_AA_standard_ASA_CWilke.txt", col_names = T)
#### files with both chains of the complex
pdb_dir <- "PyMOL_figures/pdbs/clean"
### files with only coordinates for Gsp1 - always chain A
chainA_dir <- "PyMOL_figures/pdbs/clean_chainA/"
### files with only partner coordinates, any chain but A
other_chain_dir <- "PyMOL_figures/pdbs/clean_other_chain/"
####
files_list <- list.files(path = pdb_dir)
files_list <- files_list[grepl('pdb', files_list)]
all_chains_file_paths <- file.path(pdb_dir, files_list)
chainA_file_paths <- file.path(chainA_dir, files_list)
other_chain_file_paths <- file.path(other_chain_dir, files_list)
index <- read_tsv("Scripts/complex_structure_analyses/index.txt", col_names = T)
### make the empty tibble that will contain changes in ASA upon complex formation for every residue on every chain in the complex
interface_SASA_table <- tibble(file_path = character(), chain = character(), resno = integer(), resid = character(),
                               interface = character(), rosetta_num = integer(), deltarASA = double(), deltaASA = double(), rASAc = double(),
                               aa = character())
### call dssp from bio3d and use definitions from Levy E, JMB, 2010 for interface core, support and rim
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
    mutate('deltaASA' = SASA_A_monomer - SASA_A_complex) %>% 
    mutate("interface" = ifelse((rASAm >= 0.25 & rASAc < 0.25 & deltarASA > 0), "core", "rim")) %>% 
    mutate("interface" = ifelse(rASAm < 0.25, "support", interface)) %>% 
    mutate('interface' = ifelse(deltarASA == 0, 'not_interface', interface)) %>% 
    select(file_path, chain, resno, resid, deltarASA, deltaASA, rASAc, interface) %>% 
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
    mutate('deltaASA' = SASA_B_monomer - SASA_B_complex) %>% 
    mutate("interface" = ifelse((rASAm >= 0.25 & rASAc < 0.25 & deltarASA > 0), "core", "rim")) %>% 
    mutate("interface" = ifelse(rASAm < 0.25, "support", interface)) %>% 
    mutate('interface' = ifelse(deltarASA == 0, 'not_interface', interface)) %>% 
    select(file_path, chain, resno, resid, deltarASA, deltaASA, rASAc, interface) %>% 
    bind_cols(., unique(tibble('aa' = pdbseq(pdbB), 'pdb_seq_number' = names(pdbseq(pdbB))))) %>% 
    filter(resno > 0)
  length <- length(chainB_deltaSASA$aa)
  chainB_deltaSASA <- chainB_deltaSASA %>% mutate('rosetta_num' = seq(1, length, 1))
  interface_SASA_table <- interface_SASA_table %>% 
    bind_rows(., select(chainA_deltaSASA, -pdb_seq_number)) %>% 
    bind_rows(., select(chainB_deltaSASA, -pdb_seq_number))
}

#### combine the tibble that contains SASA and interface information with the index table
interface_SASA_table <- index %>%
  select(file_path, chain, protein, partner, species) %>% 
  inner_join(., interface_SASA_table, by = c('file_path', 'chain'))
files <- interface_SASA_table %>% pull(file_path) %>% unique()

### remove NUP60, it's the same thing as NUP1
interface_SASA_table <- interface_SASA_table %>% 
  filter(! partner == 'NUP60')
#### index file contained the sequences of yeast homologs for each of the proteins (whether the actual structure is from yeast or not)
### use seqbind and seqaln from bio3d to align the pdb sequence with the yeast sequence
alignments_with_yeast_tibble <- tibble('file_path' = character(), 'chain' = character(), 'pdb_seq' = character(), 'yeast_seq' = character()) 
for (i in seq_along(files)) {
  file_path_temp <- files[i]
  temp <- interface_SASA_table %>% filter(file_path == file_path_temp)
  temp_chains <- temp %>% pull(chain) %>% unique()
  for (j in seq_along(temp_chains)) {
    temp_chain <- temp_chains[j]
    seq <- temp %>% 
      filter(chain == temp_chain) %>%
      pull(aa) %>%
      str_c(collapse = '')
    index_temp <- index %>% 
      filter(file_path == file_path_temp & chain == temp_chain)
    sequences <- seqbind(seq, index_temp$protein_sequence[1])
    full_seq_aln <- seqaln(sequences) 
    aln_tibble <- tibble('file_path' = file_path_temp, 'chain' = temp_chain, 'pdb_seq' = full_seq_aln$ali[1,], 
                         'yeast_seq' = full_seq_aln$ali[2,]) %>% 
                  add_column(yeast_num = 0, rosetta_num = 0)
      
    yeast_counter <- 0
    rosetta_counter <- 0
    for (k in seq(1, length(aln_tibble$file_path))) {
      if (aln_tibble[k,]$yeast_seq != '-') {
        yeast_counter <- yeast_counter + 1
        aln_tibble[k,]$yeast_num <- yeast_counter
        }
      if (aln_tibble[k,]$pdb_seq != '-') {
        rosetta_counter <- rosetta_counter + 1
        aln_tibble[k,]$rosetta_num <- rosetta_counter
        }
    }
    # aln_tibble <- tibble('file_path' = file_path_temp, 'chain' = temp_chain, 'pdb_seq' = full_seq_aln$ali[1,], 'yeast_seq' = full_seq_aln$ali[2,]) %>%
    #   filter(! yeast_seq == '-') %>% 
    #   mutate('yeast_resnum' = seq(1, length(full_seq_aln$ali[2,]), 1)) %>% 
    #   filter(! pdb_seq == '-')
    # length <- length(aln_tibble$pdb_seq)
    # aln_tibble <- aln_tibble %>% 
    #   mutate('rosetta_num' = seq(1, length, 1))
    alignments_with_yeast_tibble <- alignments_with_yeast_tibble %>% 
      bind_rows(., aln_tibble)
  }
}

#### merge and export the file
merged <- interface_SASA_table %>% 
  inner_join(., alignments_with_yeast_tibble, by = c('file_path', 'chain', 'rosetta_num', 'aa' = 'pdb_seq')) %>% 
  mutate('identical' = ifelse(aa == yeast_seq, TRUE, FALSE))

write_tsv(merged, path = 'Data/Gsp1_interfaces_SASA_and_conservation.txt')


