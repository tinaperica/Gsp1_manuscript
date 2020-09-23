run ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/set_ucsf_colors_in_pymol.pml

bg_color white

load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/3m1i.pdb_flsq.pdb"


create RanGTP, 3m1i.pdb_flsq and chain A and resi 1-182
create GTP, 3m1i.pdb_flsq and chain A and not polymer
delete 3m1i.pdb_flsq

as cartoon
color ucsf_navy2, RanGTP

util.cbaw not polymer
show sticks, not polymer

select RanGTP_mutations, RanGTP and resi 34+58+78+79+80+84+101+102+105+108+112+115+129+132+137+139+141+143+147+148+154+157+169+180 and name ca
show sphere, RanGTP_mutations

deselect

select RanGTP_Switch1, RanGTP and resi 41-47
select RanGTP_Switch2, RanGTP and resi 69-77
color ucsf_blue3, RanGTP_Switch1
color ucsf_pink3, RanGTP_Switch2

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


png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Revisions/Main Figures/Figure1/Fig1B_mutations_and_switch_front.png, dpi = 1200, 0, 0, -1

rotate y, 180

png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Revisions/Main Figures/Figure1/Fig1B_mutations_and_switch_back.png, dpi = 1200, 0, 0, -1
