set_color ucsf_navy1, [0.02, 0.13, 0.29]
set_color ucsf_navy2, [0.31, 0.39, 0.5]
set_color ucsf_navy3, [0.61, 0.65, 0.71]
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

##### load Ran, YRB1
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/3m1i.pdb_flsq.pdb"
create Ran, 3m1i.pdb_flsq and chain A and not res 182-1178
create GTP, 3m1i.pdb_flsq and chain A and res 1177-1178
create YRB1, 3m1i.pdb_flsq and chain B
create CRM1, 3m1i.pdb_flsq and chain C
delete 3m1i.pdb_flsq
color ucsf_navy2, Ran
util.cbaw GTP
color ucsf_yellow1, YRB1

##### load CRM1
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/3m1i1.pdb_flsq.pdb"
create CRM1, 3m1i1.pdb_flsq and chain C
delete 3m1i1.pdb_flsq
color ucsf_green1, CRM1

##### load RanGEF
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1i2m.pdb_flsq.pdb"
create RanGEF, 1i2m.pdb_flsq and chain B
delete 1i2m.pdb_flsq
color ucsf_cyan1, RanGEF

##### load RanGAP
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1k5d1.pdb_flsq.pdb"
create RanGAP, 1k5d1.pdb_flsq and chain C
delete 1k5d1.pdb_flsq
color ucsf_orange1, RanGAP

##### load KAP95
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/2bku.pdb_flsq.pdb"
create KAP95, 2bku.pdb_flsq and chain B
delete 2bku.pdb_flsq
color ucsf_green1, KAP95

##### load LOS1
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/3icq.pdb_flsq.pdb"
create LOS1, 3icq.pdb_flsq and chain T
delete 3icq.pdb_flsq
color ucsf_green1, LOS1

##### load MTR10
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/4ol0.pdb_flsq.pdb"
create MTR10, 4ol0.pdb_flsq and chain B
delete 4ol0.pdb_flsq
color ucsf_green1, MTR10

##### load CSE1
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1wa51.pdb_flsq.pdb"
create CSE1, 1wa51.pdb_flsq and chain C
delete 1wa51.pdb_flsq
color ucsf_green1, CSE1

##### load MSN5
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/3a6p.pdb_flsq.pdb"
create MSN5, 3a6p.pdb_flsq and chain F
delete 3a6p.pdb_flsq
color ucsf_green1, MSN5

##### load SRP1
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1wa5.pdb_flsq.pdb"
create SRP1, 1wa5.pdb_flsq and chain B
delete 1wa5.pdb_flsq
color ucsf_pink2, SRP1

##### load PSE1
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/3w3z.pdb_flsq.pdb"
create PSE1, 3w3z.pdb_flsq and chain C
delete 3w3z.pdb_flsq
color ucsf_green1, PSE1

##### load NTF2
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1a2k.pdb_flsq.pdb"
create NTF2, 1a2k.pdb_flsq and chain B
delete 1a2k.pdb_flsq
color ucsf_purple1, NTF2

##### load KAP104
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1qbk.pdb_flsq.pdb"
create KAP104, 1qbk.pdb_flsq and chain B
delete 1qbk.pdb_flsq
color ucsf_green1, KAP104

##### load YRB2
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/3wyf.pdb_flsq.pdb"
create YRB2, 3wyf.pdb_flsq and chain B
delete 3wyf.pdb_flsq
color ucsf_yellow1, YRB2

##### Set view
set_view (\
   -0.310891718,    0.697577894,   -0.645549059,\
   -0.647396207,   -0.652700663,   -0.393524557,\
   -0.695863664,    0.295580894,    0.654529870,\
    0.000000000,    0.000000000, -500.342407227,\
   55.409477234,  -18.293485641,   70.340736389,\
  393.955810547,  606.729003906,  -20.000000000 )
hide all
show surface, Ran
show stick, GTP
bg_color white
ray
set opaque_background, 0
zoom vis

###### Print uncolored Ran as cartton
hide surface, Ran
show cartoon, Ran
color ucsf_navy2, Ran
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/Ran_cartoon_for_globe.png, dpi = 1200, 0, 0, -1
hide cartoon, Ran
show surface, Ran
color ucsf_navy2, Ran


###### Print uncolored
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/uncolored_front.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/uncolored_back.png, dpi = 1200, 0, 0, -1
rotate y, 180  


##### Color based on YRB1 or YRB1 overlap
select br. Ran and  ((all within 3 of YRB1) or (all within 3 of YRB2))
color ucsf_yellow1, sele
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/YRB1_front.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/YRB1_back.png, dpi = 1200, 0, 0, -1
rotate y, 180  
color ucsf_navy2, sele

##### Color based on karyopherin overlap
select br. Ran and ((all within 3 of CRM1) or (all within 3 of KAP95) or (all within 3 of LOS1) or (all within 3 of MTR10) or (all within 3 of CSE1) or (all within 3 of MSN5) or (all within 3 of PSE1) or (all within 3 of KAP104))
color ucsf_green1, sele
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/karyopherins_front.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/karyopherins_back.png, dpi = 1200, 0, 0, -1
rotate y, 180  
color ucsf_navy2, sele

##### Color based on RanGEF overlap
select br. Ran and all within 3 of RanGEF
color ucsf_cyan1, sele
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/GEF_front.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/GEF_back.png, dpi = 1200, 0, 0, -1
rotate y, 180  
color ucsf_navy2, sele

##### Color based on RanGAP overlap
select br. Ran and all within 3 of RanGAP
color ucsf_orange1, sele
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/GAP_front.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/GAP_back.png, dpi = 1200, 0, 0, -1
rotate y, 180  
color ucsf_navy2, sele

##### Color based on SRP1 overlap
select br. Ran and all within 3 of SRP1
color ucsf_pink2, sele
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/SRP1_front.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/SRP1_back.png, dpi = 1200, 0, 0, -1
rotate y, 180  
color ucsf_navy2, sele


##### Color based on NTF2 overlap
select br. Ran and all within 3 of NTF2
color ucsf_purple1, sele
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/NTF2_front.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/interfaces_colored/NTF2_back.png, dpi = 1200, 0, 0, -1
rotate y, 180  
color ucsf_navy2, sele



