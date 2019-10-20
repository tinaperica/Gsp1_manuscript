run "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/set_ucsf_colors_in_pymol.pml"

bg_color white
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/3m1i.pdb_flsq.pdb"
create Ran, chain A and resi 1-181+1177
delete 3m1i.pdb_flsq
select mutations, Ran and resi 112
color ucsf_navy2, Ran
as surface
set transparency, 0.15
color ucsf_pink1, mutations
show spheres, mutations
color black, not polymer
show sticks, not polymer
set_view (\
   -0.310891718,    0.697577894,   -0.645549059,\
   -0.647396207,   -0.652700663,   -0.393524557,\
   -0.695863664,    0.295580894,    0.654529870,\
    0.000000000,    0.000000000, -500.342407227,\
   55.409477234,  -18.293485641,   70.340736389,\
  393.955810547,  606.729003906,  -20.000000000 )
bg_color white
ray
set opaque_background, 0
zoom vis
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/show_R112.png, dpi = 1200, 0, 0, -1
