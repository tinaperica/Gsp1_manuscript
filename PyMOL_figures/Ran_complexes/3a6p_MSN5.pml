run "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/set_ucsf_colors_in_pymol.pml"

bg_color white
load "/Users/cjmathy/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/3a6p.pdb_flsq.pdb"
create Ran, chain A
create MSN5, chain F
delete 3a6p.pdb_flsq
color ucsf_navy1, Ran  
color ucsf_green1, MSN5
color gray70, not polymer
as cartoon
show sticks, not polymer
set_view (\
   -0.310891718,    0.697577894,   -0.645549059,\
   -0.647396207,   -0.652700663,   -0.393524557,\
   -0.695863664,    0.295580894,    0.654529870,\
    0.000000000,    0.000000000, -500.342407227,\
   55.409477234,  -18.293485641,   70.340736389,\
  393.955810547,  606.729003906,  -20.000000000 )
png /Users/cjmathy/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/MSN5.png, dpi = 1200, 0, 0, -1
