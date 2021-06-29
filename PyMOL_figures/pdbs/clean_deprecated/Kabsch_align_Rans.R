#read in all the individual Ran structures
library(tidyverse)
library(bio3d)
files <- read_tsv("PyMOL_figures/pdbs/structures.lst", col_names = F)$X1
files <- files[grepl('pdb', files)]
# read all the pdb files in one list object
pdbs <- pdbaln(files)
core <- core.find(pdbs)
plot(core)
# gives indices of residues that will be aligned, like in sievefit
core.inds <- print(core, vol = 0.5)
# sievefit all the structures and output them to the dir defined in outpath (makes the dir also)
super.rans <- pdbfit(pdbs, core.inds, outpath = "PyMOL_figures/pdbs/corefit_structures_all")


