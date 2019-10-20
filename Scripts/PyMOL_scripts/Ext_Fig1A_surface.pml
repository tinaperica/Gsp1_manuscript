run ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/set_ucsf_colors_in_pymol.pml

bg_color white

load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1ibr.pdb_flsq.pdb"


create RanGTP, 1ibr.pdb_flsq and chain A and resi 1-181
create GTP, 1ibr.pdb_flsq and chain A and not polymer
delete 1ibr.pdb_flsq

as surface
color ucsf_navy2, RanGTP
util.cbaw not polymer
show sticks, not polymer


set_view (\
    -0.307377636,    0.700903058,   -0.643627584,\
    -0.645513952,   -0.650535166,   -0.400146306,\
    -0.699165285,    0.292473137,    0.652403951,\
    -0.000099406,    0.000011697, -149.729415894,\
    54.941619873,   -3.502172470,   78.810142517,\
    43.367378235,  256.140472412,  -20.000000000 )

png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Extended_Figures/Ext_Fig1A_surface.png, dpi = 1200, 0, 0, -1