run ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/set_ucsf_colors_in_pymol.pml
bg_color white
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1i2m.pdb_flsq.pdb"
create RanGEF, chain B
delete 1i2m.pdb_flsq
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1k5d1.pdb_flsq.pdb"
create RanGAP, 1k5d1.pdb_flsq and chain C
create Ran, 1k5d1.pdb_flsq and chain A and resi 1-182+1250-1251
delete 1k5d1.pdb_flsq
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1k5d.pdb_flsq.pdb"
create YRB1, 1k5d.pdb_flsq and chain B
delete 1k5d.pdb_flsq
color ucsf_navy1, Ran  
color ucsf_cyan1, RanGEF
color ucsf_orange1, RanGAP
color ucsf_yellow1, YRB1
# color ucsf_pink1, not polymer
util.cbay resn GNP
set valence, off
as cartoon
show sticks, not polymer
select T34, Ran and resi 32
select GAP_int, Ran and resi 130
select GEF_int, Ran and resi 99+103+106
select all_mutations, Ran and resi 32+56+76+77+78+82+99+100+103+106+110+113+127+130+135+137+139+141+145+146+152+155+167+178
select apms_mutations, Ran and resi 32+56+76+77+78+99+103+106+110+130+139+141+145+146+155+178
show spheres, apms_mutations
color ucsf_green2, apms_mutations
color ucsf_cyan1, GEF_int
color ucsf_orange1, GAP_int
color ucsf_pink2, T34
set_view (\
    -0.857339561,    0.190195411,    0.478322387,\
    -0.486194849,    0.005983099,   -0.873831570,\
    -0.169060990,   -0.981727660,    0.087342508,\
     0.000000000,    0.000000000, -313.011199951,\
    45.412685394,  -10.967408180,   78.646949768,\
   246.780426025,  379.241973877,  -20.000000000 )
hide everything, RanGEF or YRB1
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/Fig2_T34_RanGAP_2.png, dpi = 1200, 0, 0, -1
hide everything, RanGAP or YRB1
show cartoon, RanGEF
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/Fig2_T34_RanGEF_2.png, dpi = 1200, 0, 0, -1
hide everything, RanGEF or RanGAP
show cartoon, YRB1
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/Fig2_T34_YRB1_2.png, dpi = 1200, 0, 0, -1
show cartoon YRB1 or RanGEF or RanGAP
rotate x, 180
hide everything, RanGAP or YRB1
show cartoon, RanGEF
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/Fig2_T34_RanGEF_2_180.png, dpi = 1200, 0, 0, -1
hide everything, RanGEF or RanGAP
show cartoon, YRB1
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/Fig2_T34_YRB1_2_180.png, dpi = 1200, 0, 0, -1
