run ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/set_ucsf_colors_in_pymol.pml

bg_color white
### all of the structures used are mammalian, so OK to use mammalian residue numbering everywhere
### that is the reason I use 1k5d for Yrb1 (although we normally use 3m1i for other analysis - 3m1i is a yeast structure)
### GEF (Srm1)
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1i2m.pdb_flsq.pdb"
create RanSrm1, 1i2m.pdb_flsq and chain A and resi 1-181
### GAP (Rna1)
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1k5d.pdb_flsq.pdb"
create RanRna1, 1k5d.pdb_flsq and chain A and resi 1-181
create Rna1nucleotide, 1k5d.pdb_flsq and not polymer
### same structure for Yrb1 and for Rna1
delete 1k5d.pdb_flsq
# Ntf2
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1a2k.pdb_flsq.pdb"
create RanNtf2, 1a2k.pdb_flsq and chain A and resi 1-181
create Ntf2nucleotide, 1a2k.pdb_flsq and not polymer
delete 1a2k.pdb_flsq
# Srp1
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1wa5.pdb_flsq.pdb"
create RanSrp1, 1wa5.pdb_flsq and chain A and resi 1-181
create Srp1nucleotide, 1wa5.pdb_flsq and not polymer
delete 1wa5.pdb_flsq
### karyopherins
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/2bku.pdb_flsq.pdb"
create Rankaryopherins, 2bku.pdb_flsq and chain A and resi 1-181
create karyopherinnucleotide, 2bku.pdb_flsq and not polymer
delete 2bku.pdb_flsq

### make a selection for all of them
select allRans, RanSrm1 or RanRna1 or RanNtf2 or RanSrp1 or Rankaryopherins
as cartoon
color ucsf_navy2, allRans
select nucleotide, not polymer
util.cbaw not polymer
select mutations, allRans and resi 32+56+76+77+78+82+99+100+103+106+110+113+127+130+135+137+139+141+145+146+152+155+167+178 and name ca

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
hide all

select GEF_core_1, RanSrm1 and resi 19+67+71+72+73+94+95+96+99+102+103+106+107+137
select GEF_other_1, RanSrm1 and resi 18+20+23+24+68+69+70+74+76+77+91+93+97+98+100+110+130+133+134+138+140+141
show cartoon, RanSrm1
show spheres, mutations and RanSrm1
show sticks, not polymer and RanSrm1
color ucsf_cyan1, GEF_core_1
color ucsf_cyan3, GEF_other_1
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/GEF.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/GEF_180.png, dpi = 1200, 0, 0, -1
rotate y, 180
hide all
select GAP_core_1, RanRna1 and resi 19+20+39+41+43+45+69+70+95+96+128+129+130
select GAP_other_1, RanRna1 and resi 18+38+40+44+46+67+68+71+72+74+75+76+79+91+93+94+97+99+100+123+126+127+132+135
show cartoon, RanRna1
show sticks, Rna1nucleotide
show spheres, mutations and RanRna1
show sticks, not polymer and RanRna1
color ucsf_orange1, GAP_core_1
color ucsf_orange3, GAP_other_1
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/GAP.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/GAP_180.png, dpi = 1200, 0, 0, -1
rotate y, 180
hide all
color ucsf_navy2, allRans
select Ntf2_core_1, RanNtf2 and resi 19+43+67+71+72+73+76+78
select Ntf2_other_1, RanNtf2 and resi 39+40+41+42+44+68+69+70+74+77
show cartoon, RanNtf2
show sticks, Ntf2nucleotide
show spheres, mutations and RanNtf2
show sticks, not polymer and RanNtf2
color ucsf_purple1, Ntf2_core_1
color ucsf_purple3, Ntf2_other_1
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/NTF2.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/NTF2_180.png, dpi = 1200, 0, 0, -1
rotate y, 180
hide all
select Yrb12_core_1, RanRna1 and resi 31+32+51+53+154+157+171+177+178+201+210+211+212
select Yrb12_core_2, RanRna1 and resi 33+34+35+55+56+57+158+169+176+179+180+181+184+204+205+209
color ucsf_yellow1, Yrb12_core_2
color ucsf_yellow3, Yrb12_core_1
show cartoon, RanRna1
show sticks, Rna1nucleotide
show spheres, mutations and RanRna1
show sticks, not polymer and RanRna1
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/YRB.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/YRB_180.png, dpi = 1200, 0, 0, -1
rotate y, 180
hide all
select Srp1_core_1, RanSrp1 and resi 95+99+130+135
select Srp1_other_1, RanSrp1 and resi 96+103+127+128+129+132+134+137
show cartoon, RanSrp1
show sticks, Srp1nucleotide
show spheres, mutations and RanSrp1
show sticks, not polymer and RanSrp1
color ucsf_pink1, Srp1_core_1
color ucsf_pink3, Srp1_other_1
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/SRP1.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/SRP1_180.png, dpi = 1200, 0, 0, -1
rotate y, 180
hide all
select karyopherins_core_1, Rankaryopherins and resi 19+20+33+44+46+70+72+73+93+96+100+123+129+130+132+135+137+143
select karyopherins_core_2, Rankaryopherins and resi 102+103+125+126+139+141+152+153+158
select karyopherins_core_3, Rankaryopherins and resi 35+64+94+107+133+154+155+166
select karyopherins_core_4, Rankaryopherins and resi 76+106+140+145
select karyopherins_core_5, Rankaryopherins and resi 77+159
select karyopherins_core_6, Rankaryopherins and resi 47+111
select karyopherins_core_7, Rankaryopherins and resi 74+75+78+81+110
show cartoon, Rankaryopherins
show sticks, karyopherinnucleotide
show spheres, mutations and Rankaryopherins
show sticks, not polymer and Rankaryopherins
color ucsf_green3, karyopherins_core_1
color ucsf_green2, karyopherins_core_2 or karyopherins_core_3
color ucsf_green1, karyopherins_core_4 or karyopherins_core_5 or karyopherins_core_6 or karyopherins_core_7
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/karyopherins.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Pymol_figures/Ran_interface_patches/karyopherins_180.png, dpi = 1200, 0, 0, -1
rotate y, 180

