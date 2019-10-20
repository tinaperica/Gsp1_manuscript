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
create Ran, chain A and resi 1-181+1177-1178
delete 3m1i.pdb_flsq
select mutations, Ran and resi 34+58+78+79+80+84+101+102+105+108+112+115+129+132+137+139+141+143+147+148+154+157+169+180
color ucsf_navy2, Ran
as cartoon
select P_loop, Ran and resi 17-23
color ucsf_yellow1, P_loop
select switch1, Ran and resi 39-45
select switch2, Ran and resi 67-75
color ucsf_green1, switch1
color ucsf_cyan1, switch2
select nucleotide, not polymer
color ucsf_pink1, nucleotide
show sticks, nucleotide
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
png ~/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_switch/show_Ran_switches.png, dpi = 1200, 0, 0, -1
