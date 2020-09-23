run ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/set_ucsf_colors_in_pymol.pml
bg_color white
set transparency, 0.6
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1i2m.pdb_flsq.pdb"
create RanGEF, chain B
create Gsp1GEF, chain A and resi 8-177
delete 1i2m.pdb_flsq
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1k5d1.pdb_flsq.pdb"
create RanGAP, 1k5d1.pdb_flsq and chain C
create Gsp1GAP, 1k5d1.pdb_flsq and chain A and resi 1-182+1250-1251
delete 1k5d1.pdb_flsq
#load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1k5d.pdb_flsq.pdb"
#create YRB1, 1k5d.pdb_flsq and chain B
#delete 1k5d.pdb_flsq
#set_view (\
 #   -0.307377636,    0.700903058,   -0.643627584,\
  #  -0.645513952,   -0.650535166,   -0.400146306,\
   # -0.699165285,    0.292473137,    0.652403951,\
   # -0.000099406,    0.000011697, -149.729415894,\
   # 54.941619873,   -3.502172470,   78.810142517,\
   # 43.367378235,  256.140472412,  -20.000000000 )
#set_view (\
 #  -0.680408180,    0.384683847,    0.623747587,\
 #  -0.516111910,   -0.855797648,   -0.035197612,\
 #  0.520261705,   -0.345870286,    0.780832946,\
 #  0.000000000,    0.000000000, -313.011199951,\
 #  45.412685394,  -10.967408180,   78.646949768,\
 #  246.780426025,  379.241973877,  -20.000000000 )
set_view (\
    -0.848093987,   -0.337605745,    0.408352822,\
    -0.443224549,    0.874349892,   -0.197646782,\
    -0.290314972,   -0.348618954,   -0.891167343,\
    -0.000000000,    0.000000000, -221.298660278,\
    47.925292969,   -9.041076660,   78.291770935,\
   174.473541260,  268.123779297,  -20.000000000 )
select T34, Gsp1GAP and resi 32
select GEF_core_1, (Gsp1GEF) and resi 19+67+71+72+73+94+95+96+99+102+103+106+107+137 and not name c+n+o
#select GEF_other_1, (Gsp1GEF) and resi 18+20+23+24+68+69+70+74+76+77+91+93+97+98+100+110+130+133+134+138+140+141 and name ca
select GAP_core_1, (Gsp1GAP) and resi 19+20+39+41+43+45+69+70+95+96+128+129+130 and not name c+n+o
#select GAP_other_1, (Gsp1GAP) and resi 18+38+40+44+46+67+68+71+72+74+75+76+79+91+93+94+97+99+100+123+126+127+132+135 and name ca
#select all_mutations, (Gsp1GEF or Gsp1GAP) and resi 32+56+76+77+78+82+99+100+103+106+110+113+127+130+135+137+139+141+145+146+152+155+167+178 and name ca
select apms_mutations, (Gsp1GEF or Gsp1GAP) and resi 32+56+76+77+78+99+103+106+110+130+139+141+145+146+155+178
#select GAP_int, GAP_core_1 or GAP_other_1
#select GEF_int, GEF_core_1 or GEF_other_1
color ucsf_navy2, Gsp1GEF
color ucsf_navy2, Gsp1GAP
color ucsf_cyan1, RanGEF
color ucsf_orange1, RanGAP
#color ucsf_yellow1, YRB1
util.cbay resn GNP
set valence, off
hide all
show cartoon, Gsp1GEF
#show cartoon, Gsp1GAP
show spheres, GEF_core_1 and apms_mutations
color ucsf_cyan1, GEF_core_1 and apms_mutations
show surface, RanGEF
color ucsf_pink2, T34
set sphere_transparency=0.5, T34
show spheres, T34
util.cbaw not polymer
show sticks, not polymer
set valence, off
set depth_cue, off
set specular, off
set ray_shadows, off
set surface_quality, 3
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Revisions/Main Figures/Figure2/Fig2_T34_RanGEF_2.0.png, dpi = 1200, 0, 0, -1

hide all
#color ucsf_navy2, Gsp1GEF
color ucsf_navy2, Gsp1GAP
#show cartoon, Gsp1GEF
show cartoon, Gsp1GAP
show surface, RanGAP
show spheres, GAP_core_1 and apms_mutations
color ucsf_orange1, GAP_core_1 and apms_mutations
color ucsf_pink2, T34
show spheres, T34
util.cbaw not polymer
show sticks, not polymer
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Revisions/Main Figures/Figure2/Fig2_T34_RanGAP_2.0.png, dpi = 1200, 0, 0, -1

#hide all
#show cartoon, Gsp1GEF
#show cartoon, Gsp1GAP
#show spheres, GAP_int
#color ucsf_orange3, GAP_int
#color ucsf_orange1, GAP_int and apms_mutations
#color ucsf_pink2, T34
#show spheres, T34
#util.cbaw not polymer
#show sticks, not polymer

#png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Revisions/Main Figures/Figure2/Fig2_T34_RanGAP.png, dpi = 1200, 0, 0, -1


#hide everything, RanGEF or YRB1
#png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/Fig2_T34_RanGAP_2.png, dpi = 1200, 0, 0, -1
#hide everything, RanGAP or YRB1
#show cartoon, RanGEF
#png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/Fig2_T34_RanGEF_2.png, dpi = 1200, 0, 0, -1
#hide everything, RanGEF or RanGAP
#show cartoon, YRB1
#png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/Fig2_T34_YRB1_2.png, dpi = 1200, 0, 0, -1
#show cartoon YRB1 or RanGEF or RanGAP
#rotate x, 180
#hide everything, RanGAP or YRB1
#show cartoon, RanGEF
#png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/Fig2_T34_RanGEF_2_180.png, dpi = 1200, 0, 0, -1
#hide everything, RanGEF or RanGAP
#show cartoon, YRB1
#png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_complexes/Fig2_T34_YRB1_2_180.png, dpi = 1200, 0, 0, -1
