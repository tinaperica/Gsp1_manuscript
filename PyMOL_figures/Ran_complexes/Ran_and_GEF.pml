run "~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/set_ucsf_colors_in_pymol.pml"

load "/Users/tperica/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/pdbs/clean/corefit_structures_all/1i2m.pdb_flsq.pdb"
bg_color white
create Ran, chain A
create RanGEF, chain B
delete 1i2m
color navy1, Ran  
color cyan1, GEF
select P_loop, Ran and resi 17-23
color red, P_loop
select switch1, Ran and resi 39-45
select switch2, Ran and resi 67-75
color limon, switch1
color yelloworange, switch2
select R108, Ran and resi 106
select K101, Ran and resi 99
select R78, Ran and resi 76
select D79, Ran and resi 77
select R112, Ran and resi 110
color gray70, not polymer
as cartoon
show sticks, not polymer
show sticks, R108 or K101 or R78 or D79 or R112
color hotpink, R108 or K101 or R78 or D79 or R112
select mutations, Ran and resi 32+56+76+77+78+82+99+100+103+106+110+113+127+130+135+137+139+141+145+146+152+155+167+178
