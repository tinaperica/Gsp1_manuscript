run ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/set_ucsf_colors_in_pymol.pml
bg_color white
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ras/corefit_structures_all/4m1w.pdb_flsq.pdb"
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ras/corefit_structures_all/6gj8.pdb_flsq.pdb"
load "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ras/corefit_structures_all/6oim.pdb1_flsq.pdb"
create Ras, 4m1w.pdb_flsq or 6gj8.pdb_flsq or 6oim.pdb1_flsq and not hetatm
create RasKessler, 6gj8.pdb_flsq and not hetatm
create nucleotide, 6gj8.pdb_flsq and resi 203
create Ostrem, 4m1w.pdb_flsq and resi 201
create Kessler, 6gj8.pdb_flsq and resi 201
create Canon, 6oim.pdb1_flsq and resi 303
#select alpha3_Grant_region, Ras and resi 99-108
select switchI, Ras and resi 32-38
select switchII, Ras and resi 59-69
delete 4m1w.pdb_flsq.pdb
delete 6gj8.pdb_flsq.pdb
delete 6oim.pdb1_flsq.pdb
select Gsp1T34, RasKessler and resi 25
select Gsp1R78, RasKessler and resi 68
select Gsp1D79, RasKessler and resi 69
select Gsp1H141, RasKessler and resi 134
select Gsp1Q147, RasKessler and resi 140
select Gsp1Y157, RasKessler and resi 150
hide all
util.cbaw not polymer

color ucsf_navy3, Ras  
color ucsf_green1, Kessler
color ucsf_cyan1, Ostrem
color ucsf_green3, Canon
#color pink, alpha3_Grant_region
select allosteric_in_Gsp1, Gsp1T34 or Gsp1H141 or Gsp1Q147 or Gsp1Y157
color ucsf_blue1, allosteric_in_Gsp1
color ucsf_blue3, switchI
color ucsf_pink3, switchII
show cartoon, Ras

show sticks, Ostrem or Kessler or Canon or nucleotide
show spheres, allosteric_in_Gsp1
set_view (\
     0.737445533,    0.433104068,   -0.518227816,\
    -0.669943690,    0.566275656,   -0.480059266,\
     0.085559860,    0.701225519,    0.707764566,\
     0.000000000,    0.000000000, -142.096755981,\
    13.463575363,   13.717506409,   15.233680725,\
   112.030166626,  172.163345337,  -20.000000000 )
   set valence, off
set depth_cue, off
set specular, off
set ray_shadows, off
set surface_quality, 3

png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Revisions/Extended_Figures/EDF_7/Ras_allosteric_sites.png, dpi = 1200, 0, 0, -1
rotate y, 180
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/Revisions/Extended_FIgures/EDF_7/Ras_allosteric_sites_180deg.png, dpi = 1200, 0, 0, -1
