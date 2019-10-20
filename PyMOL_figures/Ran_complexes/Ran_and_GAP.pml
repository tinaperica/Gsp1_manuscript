run "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/set_ucsf_colors_in_pymol.pml"

load "/Users/tperica/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1k5d.pdb_flsq.pdb"
bg_color white
create Ran, chain A
create YRB1, chain B
delete 1k5d.pdb_flsq
load "/Users/tperica/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1k5d1.pdb_flsq.pdb"
create RanGAP, chain C
delete 1k5d1.pdb_flsq
alter open_conf, chain = "A"
alter all, segi = ""
color ucsf_navy1, Ran  
color ucsf_orange1, RanGAP
color ucsf_yellow1, YRB1
select P_loop, Ran and resi 17-23
color ucsf_pink2, P_loop
select switch1, Ran and resi 39-45
select switch2, Ran and resi 67-75
color ucsf_green2, switch1
color ucsf_blue2, switch2
color gray70, not polymer
as cartoon
show sticks, not polymer
zoom vis
