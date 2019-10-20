load "/Users/tperica/Documents/Gsp1_bioinformatics/Ran_structures/Ran_switch/Ran.pdb"
bg_color gray80
create closed_conf, chain A
create open_conf, chain B
delete Ran
set_view (\
    -0.385181963,    0.922827721,   -0.004902451,\
    -0.473523319,   -0.193080112,    0.859357715,\
     0.792092502,    0.333330452,    0.511351287,\
     0.000000000,    0.000000000, -143.946884155,\
    13.440780640,   16.859657288,  -18.849510193,\
    15.762239456,  272.130798340,  -20.000000000 )
align closed_conf, open_conf
select tails, (open_conf and resi 180-217) or (closed_conf and resi 182-219)
remove tails
alter open_conf, chain = "A"
alter all, segi = ""
color paleyellow, open_conf
color paleyellow, closed_conf
select open_P_loop, open_conf and resi 17-23
select closed_P_loop, closed_conf and resi 19-25
color red, open_P_loop or closed_P_loop
select open_switch1, open_conf and resi 39-45
color forest, open_switch1
select closed_switch1, closed_conf and resi 41-47
color limon, closed_switch1
select open_switch2, open_conf and resi 67-75
color density, open_switch2
select closed_switch2, closed_conf and resi 69-77
color marine, closed_switch2
color violet, not polymer
as cartoon
hide closed_conf
show open_conf
show sticks, not polymer
png /Users/tperica/Documents/Gsp1_bioinformatics/Ran_structures/Ran_switch/open_conformation.png, width = 1280, height = 952, dpi = 720
show closed_conf
morph switch, open_conf, closed_conf