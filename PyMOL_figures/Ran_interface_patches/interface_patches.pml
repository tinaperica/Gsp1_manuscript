run ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/set_ucsf_colors_in_pymol.pml

bg_color white

load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1ibr.pdb_flsq.pdb"
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/3gj0.pdb_flsq.pdb"

create RanGTP, 1ibr.pdb_flsq and chain A and resi 1-181
create GTP, 1ibr.pdb_flsq and chain A and not polymer
delete 1ibr.pdb_flsq

create RanGDP, 3gj0.pdb_flsq and chain A and resi 9-181
delete 3gj0.pdb_flsq

as cartoon
color ucsf_blue1, RanGTP
color ucsf_navy3, RanGDP
util.cbaw not polymer
show sticks, not polymer


select RanGTP_mutations, RanGTP and resi 32+56+76+77+78+82+99+100+103+106+110+113+127+130+135+137+139+141+145+146+152+155+167+178 and name ca
show sphere, RanGTP_mutations

select RanGDP_mutations, RanGDP and resi 32+56+76+77+78+82+99+100+103+106+110+113+127+130+135+137+139+141+145+146+152+155+167+178 and name ca
show sphere, RanGDP_mutations


select Gsp1, RanGTP or RanGDP
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

select GEF_core_1, Gsp1 and resi 19+67+71+72+73+94+95+96+99+102+103+106+107+137
select GEF_other_1, Gsp1 and resi 18+20+23+24+68+69+70+74+76+77+91+93+97+98+100+110+130+133+134+138+140+141
color ucsf_cyan1, GEF_core_1
color ucsf_cyan3, GEF_other_1
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/GEF.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/GEF_180.png, dpi = 1200, 0, 0, -1
rotate y, 180
color ucsf_blue1, RanGTP
color ucsf_navy3, RanGDP
select GAP_core_1, Gsp1 and resi 19+20+39+41+43+45+69+70+95+96+128+129+130
select GAP_other_1, Gsp1 and resi 18+38+40+44+46+67+68+71+72+74+75+76+79+91+93+94+97+99+100+123+126+127+132+135
color ucsf_orange1, GAP_core_1
color ucsf_orange3, GAP_other_1
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/GAP.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/GAP_180.png, dpi = 1200, 0, 0, -1
rotate y, 180
color ucsf_blue1, RanGTP
color ucsf_navy3, RanGDP
select Ntf2_core_1, Gsp1 and resi 19+43+67+71+72+73+76+78
select Ntf2_other_1, Gsp1 and resi 39+40+41+42+44+68+69+70+74+77
color ucsf_purple1, Ntf2_core_1
color ucsf_purple3, Ntf2_other_1
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/NTF2.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/NTF2_180.png, dpi = 1200, 0, 0, -1
rotate y, 180
color ucsf_blue1, RanGTP
color ucsf_navy3, RanGDP
select Nup160_core_2, Gsp1 and resi 3+4+168+169
color ucsf_navy3, Nup160_core_1
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/NUP.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/NUP_180.png, dpi = 1200, 0, 0, -1
rotate y, 180
color ucsf_blue1, RanGTP
color ucsf_navy3, RanGDP
select Yrb12_core_1, Gsp1 and resi 31+32+51+53+154+157+171+177+178+201+210+211+212
select Yrb12_core_2, Gsp1 and resi 33+34+35+55+56+57+158+169+176+179+180+181+184+204+205+209
color ucsf_yellow1, Yrb12_core_2
color ucsf_yellow3, Yrb12_core_1
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/YRB.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/YRB_180.png, dpi = 1200, 0, 0, -1
rotate y, 180
color ucsf_blue1, RanGTP
color ucsf_navy3, RanGDP
select Srp1_core_1, Gsp1 and resi 95+99+130+135
select Srp1_other_1, Gsp1 and resi 96+103+127+128+129+132+134+137
color ucsf_pink1, Srp1_core_1
color ucsf_pink3, Srp1_other_1
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/SRP1.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/SRP1_180.png, dpi = 1200, 0, 0, -1
rotate y, 180
color ucsf_blue1, RanGTP
color ucsf_navy3, RanGDP
select karyopherins_core_1, Gsp1 and resi 19+20+33+44+46+70+72+73+93+96+100+123+129+130+132+135+137+143
select karyopherins_core_2, Gsp1 and resi 102+103+125+126+139+141+152+153+158
select karyopherins_core_3, Gsp1 and resi 35+64+94+107+133+154+155+166
select karyopherins_core_4, Gsp1 and resi 76+106+140+145
select karyopherins_core_5, Gsp1 and resi 77+159
select karyopherins_core_6, Gsp1 and resi 47+111
select karyopherins_core_7, Gsp1 and resi 74+75+78+81+110
color ucsf_green3, karyopherins_core_1
color ucsf_green2, karyopherins_core_2 or karyopherins_core_3
color ucsf_green1, karyopherins_core_4 or karyopherins_core_5 or karyopherins_core_6 or karyopherins_core_7
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/karyopherins.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/karyopherins_180.png, dpi = 1200, 0, 0, -1
rotate y, 180
color ucsf_blue1, RanGTP
color ucsf_navy3, RanGDP
