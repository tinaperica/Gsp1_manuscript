set_color ucsf_navy1, [0.02, 0.13, 0.29]
set_color ucsf_navy2, [0.31, 0.39, 0.5]
set_color ucsf_blue1, [0.09, 0.55, 0.80]
set_color ucsf_blue2, [0.36, 0.69, 0.86]
set_color ucsf_orange1, [0.96, 0.50, 0.14]
set_color ucsf_orange2, [0.97, 0.65, 0.40]
set_color ucsf_cyan1, [0.09, 0.64, 0.67]
set_color ucsf_cyan2, [0.36, 0.75, 0.77]
set_color ucsf_green1, [0.56, 0.74, 0.19]
set_color ucsf_green2, [0.69, 0.82, 0.44]
set_color ucsf_purple1, [0.44, 0.44, 0.70]
set_color ucsf_purple2, [0.61, 0.60, 0.79]
set_color ucsf_pink1, [0.93, 0.09, 0.28]
set_color ucsf_pink2, [0.95, 0.36, 0.5]
set_color ucsf_yellow1, [1.00, 0.86, 0.00]
set_color ucsf_yellow2, [1.00, 0.91, 0.3]
bg_color white
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/3m1i.pdb_flsq.pdb"
create Ran, chain A and resi 1-181+1177+1178
delete 3m1i.pdb_flsq
color not polymer
color ucsf_navy2, Ran
as surface
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
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/uncolored_front.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/uncolored_back.png, dpi = 1200, 0, 0, -1
rotate y, 180  
#######
####### CRM1
#######
select karyopherins, Ran and resi 35+41+72+46+48+49+76+77+78+79+80+83+98+105+108+112+113+127+134+135+142+143+144+145+147+154+155+156+157
color ucsf_green1, karyopherins
bg_color white
ray
set opaque_background, 0
zoom vis
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/karyopherins_front.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/karyopherins_back.png, dpi = 1200, 0, 0, -1
rotate y, 180
color ucsf_navy2, karyopherins
#######
####### NTF2
#######
select ntf2, Ran and resi 44+45+69+73+74+78+80
color ucsf_purple2, ntf2
bg_color white
ray
set opaque_background, 0
zoom vis
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/ntf2_front.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/ntf2_back.png, dpi = 1200, 0, 0, -1
rotate y, 180
color ucsf_navy2, ntf2
#######
####### RNA1
#######
# based on SASA core
#select rna1, Ran and resi 41+43+45+71+72+96+97+98+130+132
# by eye, based on structure
select rna1, Ran and resi 20-23+41-46+70-72+76-77+94-100+130-133
color ucsf_orange1, rna1
bg_color white
ray
set opaque_background, 0
zoom vis
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/rna1_front.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/rna1_back.png, dpi = 1200, 0, 0, -1
rotate y, 180
color ucsf_navy2, rna1
#######
####### SRM1
#######
select srm1, Ran and resi 21+73+74+75+96+97+98+101+104+105+108+139
color ucsf_cyan1, srm1
bg_color white
ray
set opaque_background, 0
zoom vis
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/srm1_front.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/srm1_back.png, dpi = 1200, 0, 0, -1
rotate y, 180
color ucsf_navy2, srm1
#######
####### SRP1
#######
select srp1, Ran and resi 97+101
color ucsf_pink1, srp1
bg_color white
ray
set opaque_background, 0
zoom vis
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/srp1_front.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/srp1_back.png, dpi = 1200, 0, 0, -1
rotate y, 180
color ucsf_navy2, srp1
#######
####### YRB1
#######
select yrb1, Ran and resi 34+35+36+57+58+59+169+170+171+180+181
color ucsf_yellow2, yrb1
bg_color white
ray
set opaque_background, 0
zoom vis
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/yrb1_front.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/yrb1_back.png, dpi = 1200, 0, 0, -1
rotate y, 180
color ucsf_navy2, yrb1
