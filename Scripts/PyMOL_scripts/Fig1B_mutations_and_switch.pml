run ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/set_ucsf_colors_in_pymol.pml

bg_color white

load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1ibr.pdb_flsq.pdb"
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/3gj0.pdb_flsq.pdb"


align 1ibr.pdb_flsq, 3gj0.pdb_flsq.pdb

create RanGTP, 1ibr.pdb_flsq and chain A and resi 1-181
create GTP, 1ibr.pdb_flsq and chain A and not polymer
delete 1ibr.pdb_flsq

create RanGDP, 3gj0.pdb_flsq and chain A and resi 9-181
# create GDP, 3gj0.pdb_flsq and chain A and not polymer
delete 3gj0.pdb_flsq

as cartoon
# color ucsf_navy1, RanGTP
color ucsf_blue1, RanGTP

# color ucsf_blue1, RanGDP
color ucsf_navy3, RanGDP
util.cbaw not polymer
show sticks, not polymer

select RanGTP_mutations, RanGTP and resi 32+56+76+77+78+82+99+100+103+106+110+113+127+130+135+137+139+141+145+146+152+155+167+178 and name ca
# color ucsf_navy1, RanGTP_mutations
color ucsf_blue1, RanGTP_mutations
show sphere, RanGTP_mutations

deselect

select RanGDP_mutations, RanGDP and resi 32+56+76+77+78+82+99+100+103+106+110+113+127+130+135+137+139+141+145+146+152+155+167+178 and name ca
# color ucsf_blue1, RanGDP_mutations
color ucsf_navy3, RanGDP_mutations
show sphere, RanGDP_mutations

deselect

select RanGTP_Switch1, RanGTP and resi 39-45
select RanGTP_Switch2, RanGTP and resi 67-75
select RanGDP_Switch1, RanGDP and resi 39-45
select RanGDP_Switch2, RanGDP and resi 67-75
color ucsf_green1, RanGTP_Switch1
color ucsf_green3, RanGDP_Switch1
color ucsf_yellow1, RanGTP_Switch2
color ucsf_yellow3, RanGDP_Switch2

deselect

set_view (\
    -0.307377636,    0.700903058,   -0.643627584,\
    -0.645513952,   -0.650535166,   -0.400146306,\
    -0.699165285,    0.292473137,    0.652403951,\
    -0.000099406,    0.000011697, -149.729415894,\
    54.941619873,   -3.502172470,   78.810142517,\
    43.367378235,  256.140472412,  -20.000000000 )

set valence, off
set depth_cue, off
set specular, off
set ray_shadows, off
set surface_quality, 3


png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Figure1_E-MAP/Plots/Fig1B_mutations_and_switch_front.png, dpi = 1200, 0, 0, -1

rotate y, 180

png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Figure1_E-MAP/Plots/Fig1B_mutations_and_switch_back.png, dpi = 1200, 0, 0, -1
