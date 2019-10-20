run "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/set_ucsf_colors_in_pymol.pml"

bg_color white
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1k5d.pdb_flsq.pdb"
create Ran, chain A
create YRB1, chain B
delete 1k5d.pdb_flsq
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1k5d1.pdb_flsq.pdb"
create RanGAP, 1k5d1.pdb_flsq and chain C
delete 1k5d1.pdb_flsq
color ucsf_navy2, Ran  
color ucsf_yellow1, YRB1
color ucsf_orange1, RanGAP
color ucsf_pink1, not polymer
select T34, Ran and resi 34
as cartoon
show sticks, not polymer
show spheres, T34
color ucsf_pink2, T34
set_view (\
   -0.310891718,    0.697577894,   -0.645549059,\
   -0.647396207,   -0.652700663,   -0.393524557,\
   -0.695863664,    0.295580894,    0.654529870,\
    0.000000000,    0.000000000, -500.342407227,\
   55.409477234,  -18.293485641,   70.340736389,\
  393.955810547,  606.729003906,  -20.000000000 )
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/GAP_YRB1_T34.png, dpi = 1200, 0, 0, -1
