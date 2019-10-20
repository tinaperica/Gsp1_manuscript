run ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/set_ucsf_colors_in_pymol.pml

bg_color white
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1k5d.pdb_flsq.pdb"
create Ran, chain A
create YRB1, chain B
delete 1k5d.pdb_flsq
color ucsf_navy1, Ran  
color ucsf_gray3, YRB1
util.cbay resn GNP
set valence, off
select T34, Ran and resi 32
as cartoon
hide cartoon, YRB1
show surface, YRB1
set transparency, 0.2, YRB1
set stick_radius, 0.5
show sticks, not polymer
show spheres, T34
color ucsf_pink2, T34
set_view (\
    -0.769523501,   -0.277580559,   -0.575140119,\
    -0.145600438,   -0.800619841,    0.581216455,\
    -0.621800959,    0.530997396,    0.575681746,\
    -0.000218101,   -0.000038806,  -95.072479248,\
    60.722293854,    1.792496681,   82.240623474,\
   -11.147430420,  201.625595093,  -20.000000000 )
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Figure3_Biophysics/Plots/3F_YRB1_T34.png, dpi = 1200, 0, 0, -1
