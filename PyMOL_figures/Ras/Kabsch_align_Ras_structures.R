# align Ras structures to show where the allosteric inhibitors are bound  
### and where the identified allosteric inhibitors are
library(tidyverse)
library(bio3d)
pdb_dir <- 'PyMOL_figures/Ras/pdbs/'
files_list <- list.files(path = pdb_dir)
file_paths <- file.path(pdb_dir, files_list)
# read all the pdb files in one list object
pdbs <- pdbaln(file_paths)
core <- core.find(pdbs)
plot(core)
# gives indices of residues that will be aligned, like in sievefit
core.inds <- print(core, vol = 0.5)
# sievefit all the structures and output them to the dir defined in outpath (makes the dir also)
super.ras <- pdbfit(pdbs, core.inds, outpath = "PyMOL_figures/Ras/corefit_structures_all")

