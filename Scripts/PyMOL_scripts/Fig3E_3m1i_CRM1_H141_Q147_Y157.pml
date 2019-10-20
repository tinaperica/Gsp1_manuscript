run ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/set_ucsf_colors_in_pymol.pml

bg_color white
load "/Users/cjmathy/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/3m1i1.pdb_flsq.pdb"
create Ran, chain A and res 1-181 or resn GTP
create CRM1, chain C
delete 3m1i1.pdb_flsq
color ucsf_navy2, Ran  
color ucsf_gray3, CRM1
util.cbay resn GTP
set valence, off
as cartoon
hide cartoon, CRM1
show surface, CRM1
set transparency, 0.2, CRM1
set stick_radius, 0.5
show sticks, not polymer
select H141, Ran and resi 141
show spheres, H141
color ucsf_pink2, H141
select Q147, Ran and resi 147
show spheres, Q147
color ucsf_pink2, Q147
select Y157, Ran and resi 157
show spheres, Y157
color ucsf_pink2, Y157
set_view (\
    -0.579937458,    0.428930640,    0.692597508,\
     0.314267993,   -0.666569471,    0.675961077,\
     0.751603544,    0.609673679,    0.251764596,\
    -0.000256464,   -0.000348836, -157.100402832,\
    51.504898071,  -10.603000641,   68.648071289,\
    37.852588654,  276.255279541,  -20.000000000 )
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Figure3_Biophysics/Plots/3E_CRM1_H141_Q147_Y157A.png, dpi = 1200, 0, 0, -1




