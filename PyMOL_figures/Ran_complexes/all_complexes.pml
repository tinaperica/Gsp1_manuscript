run "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/set_ucsf_colors_in_pymol.pml"

bg_color white
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1i2m.pdb_flsq.pdb"
set_name 1i2m.pdb_flsq, 1i2m
create RanGEF, 1i2m and chain B
color ucsf_cyan1, RanGEF
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1k5d1.pdb_flsq.pdb"
set_name 1k5d1.pdb_flsq, 1k5d1
create RanGAP, 1k5d1 and chain C
color ucsf_orange1, RanGAP
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/3m1i.pdb_flsq.pdb"
set_name 3m1i.pdb_flsq, 3m1i
create YRB1, 3m1i and chain B
color ucsf_yellow1, YRB1
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/2bku.pdb_flsq.pdb"
set_name 2bku.pdb_flsq, 2bku
create KAP95, 2bku and chain B
color ucsf_green1, KAP95
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/3icq.pdb_flsq.pdb"
set_name  3icq.pdb_flsq, 3icq
create LOS1, 3icq and chain T
color ucsf_green1, LOS1
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/4ol0.pdb_flsq.pdb"
set_name 4ol0.pdb_flsq, 4ol0
create MTR10, 4ol0 and chain B
color ucsf_green1, MTR10
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1wa51.pdb_flsq.pdb"
set_name 1wa51.pdb_flsq, 1wa51
create CSE1, 1wa51 and chain C
color ucsf_green1, CSE1
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/3m1i1.pdb_flsq.pdb"
set_name 3m1i1.pdb_flsq, 3m1i1
create CRM1, 3m1i1 and chain C
color ucsf_green1, CRM1
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/3a6p.pdb_flsq.pdb"
set_name 3a6p.pdb_flsq, 3a6p
create MSN5, 3a6p and chain F
color ucsf_green1, MSN5
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1wa5.pdb_flsq.pdb"
set_name 1wa5.pdb_flsq, 1wa5
create SRP1, 1wa5 and chain B
color ucsf_pink2, SRP1
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/3w3z.pdb_flsq.pdb"
set_name 3w3z.pdb_flsq, 3w3z
create PSE1, 3w3z and chain C
color ucsf_green1, PSE1
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1a2k.pdb_flsq.pdb"
set_name 1a2k.pdb_flsq, 1a2k
create NTF2, 1a2k and chain B
color ucsf_purple1, NTF2
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1qbk.pdb_flsq.pdb"
create KAP104, chain B
set_name 1qbk.pdb_flsq, 1qbk
create KAP104, 1qbk and chain B
color ucsf_green1, KAP104
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/3wyf.pdb_flsq.pdb"
set_name 3wyf.pdb_flsq, 3wyf
create YRB2, 3wyf and chain B
color ucsf_yellow1, YRB2
create Ran, 3m1i1 and chain A
color ucsf_navy1, Ran  
as cartoon
show sticks, not polymer
set_view (\
   -0.310891718,    0.697577894,   -0.645549059,\
   -0.647396207,   -0.652700663,   -0.393524557,\
   -0.695863664,    0.295580894,    0.654529870,\
    0.000000000,    0.000000000, -500.342407227,\
   55.409477234,  -18.293485641,   70.340736389,\
  393.955810547,  606.729003906,  -20.000000000 )

