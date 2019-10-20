load "/Users/tperica/Documents/Gsp1_bioinformatics/Ran_structures/Ran_switch/Ran.pdb"
bg_color white
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
color skyblue, open_conf
color skyblue, closed_conf
select open_P_loop, open_conf and resi 17-23
select closed_P_loop, closed_conf and resi 19-25
color red, open_P_loop or closed_P_loop
select open_switch1, open_conf and resi 39-45
select closed_switch1, closed_conf and resi 41-47
select open_switch2, open_conf and resi 67-75
select closed_switch2, closed_conf and resi 69-77
color limon, open_switch1 or closed_switch1
color yelloworange, open_switch2 or closed_switch2
select closedThr44, closed_conf and resi 44
select openThr44, open_conf and resi 42
select closedPhe37, closed_conf and resi 37
select openPhe37, open_conf and resi 35
select closedThr34, closed_conf and resi 34
select openThr34, open_conf and resi 32
select closedGln147, closed_conf and resi 147
select openGln34, open_conf and resi 145
color gray70, not polymer
as cartoon
show sticks, not polymer
select mutations_open, open_conf and resi 32+56+76+77+78+82+99+100+103+106+110+113+127+130+135+137+139+141+145+146+152+155+167+178
select mutations_closed, closed_conf and resi 34+58+78+79+80+84+101+102+105+108+112+115+129+132+137+139+141+143+147+148+154+157+169+180