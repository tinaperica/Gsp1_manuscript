library(bio3d)
library(tidyverse)
library(RColorBrewer)



standard_SASA <- read_tsv("Scripts/complex_structure_analyses/per_AA_standard_ASA_CWilke.txt", col_names = T)
switch_loops <- tibble(resi = c(19:25, 41:47, 69:77), "loop" = c(rep("ploop", 7), rep("switch1", 7), rep("switch2",9)))
pdb_dir <- "PyMOL_figures/pdbs/clean"
chainA_dir <- "PyMOL_figures/pdbs/clean_chainA"
files_list <- list.files(path = pdb_dir)
all_chains_file_paths <- file.path(pdb_dir, files_list)
chainA_file_paths <- file.path(chainA_dir, files_list)
nucleotide_contacts_table <- data.frame()
interface_SASA_table <- tibble(file_path = character(), resno = integer(), resid = character(),
              SASA_difference = double(), percent_SASA_change = double(), interface = character())

for (i in seq_along(all_chains_file_paths)) {
  file_path <- all_chains_file_paths[i]
  chainA_path <- chainA_file_paths[i]
  pdb <- read.pdb(file_path, rm.alt = FALSE)
  pdbA <- read.pdb(chainA_path, rm.alt = FALSE)
  SASA_monomer <- dssp(pdbA)$acc
  SASA_complex <- dssp(pdb)$acc[1:length(SASA_monomer)]
  coord.temp <- as_tibble(pdb$atom[, c("chain", "resno", "resid", "elety", "x", "y", "z")]) %>%
    filter(resid != "HOH" & resid != "SO4" & resid != "MG" & 
             resid != "GNP" & resid != "GDP"  & resid != "GTP" & resid != "AF3" & resid != "GOL" &
             resid != "EDO" & resid != "CL")
  chainA_deltaSASA <- coord.temp %>% 
    filter(chain == "A") %>% 
    select(resno, resid) %>% unique() %>% 
    add_column(., SASA_monomer, SASA_complex, file_path) %>%
    inner_join(., standard_SASA, by = c("resid" = "AA")) %>% 
    select(-A, -Residue) %>% 
    mutate("rASAm" = SASA_monomer/maxASA, "rASAc" = SASA_complex/maxASA) %>% 
    mutate("deltarASA" = rASAm - rASAc) %>% 
    filter(deltarASA > 0) %>% 
    mutate("interface" = ifelse((rASAm >= 0.25 & rASAc < 0.25 & deltarASA >= 0.25), "core", "rim")) %>% 
    mutate("interface" = ifelse(rASAm < 0.25, "support", interface)) %>% 
    select(file_path, resno, resid, deltarASA, interface)
  interface_SASA_table <- interface_SASA_table %>% bind_rows(., chainA_deltaSASA)
  nucleotide_coord <- as_tibble(pdb$atom[, c("chain", "resno", "resid", "elety", "x", "y", "z")]) %>%
    filter(resid != "HOH" & resid != "SO4" & resid != "AF3" & resid != "GOL" &
             resid != "EDO" & resid != "CL" & (resid == "MG" | resid == "GNP" | resid == "GDP" | resid == "GTP")) %>% 
    mutate("unique" = str_c(chain, resno, elety, sep = "_"))
  if (nrow(nucleotide_coord)) {
  nucleotide_coord_matrix <- as.matrix(nucleotide_coord[, c("x", "y", "z")])
    rownames(nucleotide_coord_matrix) <- nucleotide_coord$unique
    Gsp1_coord <- as_tibble(pdb$atom[, c("chain", "resno", "resid", "elety", "x", "y", "z")]) %>%
      filter(chain == "A" & resid != "HOH" & resid != "SO4" & resid != "MG" & 
               resid != "GNP" & resid != "GDP"  & resid != "GTP" & resid != "AF3" & resid != "GOL" &
               resid != "EDO" & resid != "CL") %>% 
      mutate("unique" = str_c(chain, resno, elety, sep = "_"))
    Gsp1_coord_matrix <- as.matrix(Gsp1_coord[, c("x", "y", "z")])
    rownames(Gsp1_coord_matrix) <- Gsp1_coord$unique
    distances <- dist.xyz(Gsp1_coord_matrix, nucleotide_coord_matrix)
    colnames(distances) <- rownames(nucleotide_coord_matrix)
    distances_df <- data.frame(distances)
    distances_df <- cbind(distances_df, "unique_Gsp1" = as.character(rownames(Gsp1_coord_matrix)))
    nucleotide_contacts <- as_tibble(distances_df) %>% 
      gather(unique_nucleotide, dist, -unique_Gsp1) %>%
      filter(dist < 5) %>%     ##### threshold for nucleotide contact is 5 A
      separate(col = unique_Gsp1, into = c("Gsp1", "resnum1", "atom1"), sep = "_", convert = T) %>%
      separate(col = unique_nucleotide, into = c("nucleotide", "resnum2", "atom2"), sep = "_", convert = T) %>%
      group_by(resnum1, resnum2) %>%
      summarise("number_of_contacts" = n()) %>%
      arrange(resnum1, resnum2)
    if (nrow(nucleotide_contacts) > 0) {
      nucleotide_contacts_table <- rbind(nucleotide_contacts_table, 
                              as_tibble(data.frame(file_path, nucleotide_contacts)))
    }
  }
  # chains <- unique(coord.temp$chain)
  # coords <- list()
  # for (ch in seq_along(chains)) {
  #   chain_id <- chains[ch]
  #   coords[[chain_id]] <- coord.temp %>% 
  #     filter(chain == chain_id) %>% 
  #     mutate("unique" = str_c(chain, resno, elety, sep = "_"))
  # }
  # 
  #pairs_of_chains <- combn(chains, 2)
  # for (p in seq_along(chains)) {
  #   chain1 <- "A"
  #   chain2 <- chains[p]
  #   if (chain2 != "A") {
  #     matrix1 <- as.matrix(coords[[chain1]][, c("x", "y", "z")])
  #     rownames(matrix1) <- coords[[chain1]]$unique
  #     matrix2 <- as.matrix(coords[[chain2]][, c("x", "y", "z")])
  #     rownames(matrix2) <- coords[[chain2]]$unique
  #     distances <- dist.xyz(matrix1, matrix2)
  #     colnames(distances) <- rownames(matrix2)
  #     distances_df <- data.frame(distances)
  #     distances_df <- cbind(distances_df, "unique_chain1" = as.character(rownames(matrix1)))
  #     contacts_df <- as_tibble(distances_df) %>% 
  #       gather(unique_chain2, dist, -unique_chain1) %>%
  #       filter(dist < 4) %>%
  #       separate(col = unique_chain1, into = c("chain1", "resnum1", "atom1"), convert = T) %>%
  #       separate(col = unique_chain2, into = c("chain2", "resnum2", "atom2"), convert = T) %>%
  #       group_by(resnum1, resnum2) %>%
  #       summarise("number_of_contacts" = n()) %>%
  #       arrange(resnum1, resnum2)
  #     if (nrow(contacts_df) > 0) {
  #       contacts_table <- rbind(contacts_table, as_tibble(data.frame(file_path, chain1, chain2, contacts_df)))
  #     }
  # }
  #}
}

#contacts_table <- as_tibble(contacts_table)
#write_delim(contacts_table, "4A_contacts.txt", delim = "\t")
index <- read_tsv("index.txt", col_names = T)
Ran_seq_index <- read_tsv("Gsp1_to_Ran.index", col_names = T)

### residues contacting the nucleotide (based on the 6 A threshold)
nuc_merged <- inner_join(nucleotide_contacts_table, index, by = c("file_path"))
nuc_merged_yeast <- nuc_merged %>%
  filter(species == "yeast") %>%
  select(partner, "yeastresnum" = resnum1, number_of_contacts, pdb_id)
nuc_merged_mamm <- nuc_merged %>%
  filter(species == "mammalian") %>%
  inner_join(., Ran_seq_index, by = c("resnum1" = "pdb_seq_num")) %>%
  select(partner, "yeastresnum" = ref_seq_num, number_of_contacts, pdb_id)
nucleotide_currated <- bind_rows(nuc_merged_yeast, nuc_merged_mamm)
write_tsv(nucleotide_currated, "Gsp1_residues_in_contact_with_nucleotides.txt")
#############################################################################


######## Gsp1 interfaces based on delta_rASA
merged_SASA <- interface_SASA_table %>% 
  inner_join(., index, by = "file_path") %>% 
  select(-chain)
merged_yeast_SASA <- merged_SASA %>% 
  filter(species == "yeast") %>% 
  select(partner, "yeastresnum" = resno, deltarASA, interface, pdb_id) %>% 
  mutate("deltarASA" = round(deltarASA, 2))
merged_mamm_SASA <- merged_SASA %>%
  filter(species == "mammalian") %>%  
  inner_join(., Ran_seq_index, by = c("resno" = "pdb_seq_num")) %>% 
  select(partner, "yeastresnum" = ref_seq_num, deltarASA, interface, pdb_id) %>% 
  mutate("deltarASA" = round(deltarASA, 2))
currated_SASA <- bind_rows(merged_yeast_SASA, merged_mamm_SASA) %>% 
  group_by(yeastresnum, partner) %>% 
  filter(deltarASA == max(deltarASA, na.rm = T)) %>% 
  filter(deltarASA > 0.05) %>% 
  select(-pdb_id) %>% 
  arrange(interface, yeastresnum, partner) %>% 
  unique() %>% 
  ungroup()
write_tsv(currated_SASA, "SASA_interfaces.txt")

