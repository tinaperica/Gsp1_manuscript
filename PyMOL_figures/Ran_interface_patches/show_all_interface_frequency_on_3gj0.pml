set_color ucsf_pink1, [0.93, 0.09, 0.28]
set_color ucsf_pink2, [0.95, 0.36, 0.50]
set_color ucsf_pink3, [0.97, 0.64, 0.72]
set_color ucsf_purple1, [0.44, 0.44, 0.70]
set_color ucsf_purple2, [0.61, 0.60, 0.79]
set_color ucsf_purple3, [0.78, 0.78, 0.88]
bg_color white
select chelix, 3gj0 and resi 190-207
remove chelix
### in zero interfaces
select zerointerfaces, 3gj0 and resi 5+6+9+10+13+14+15+16+17+22+25+26+27+28+61+63+65+66+80+83+85+86+87+88+89+90+92+101+108+112+115+116+117+118+119+120+121+122+131+149+150+151+160+164+170+187+188+189+190+191+192+193+194+195+196+199+215+216+217+218+219
color white, zerointerfaces
### in 1 interface
select oneinterface, 3gj0 and resi 7+21+23+24+36+98+105+161+165+186+197+198+200+202+213+214
color ucsf_pink3, oneinterface
### in 2-4 interfaces
select in2to4int, 3gj0 and resi 1+2+3+4+8+11+18+20+30+31+32+33+34+42+48+49+50+51+52+53+54+55+56+57+58+59+60+67+73+84+97+102+104+109+114+123+124+125+132+133+135+136+138+144+148+157+158+162+167+168+169+171+172+174+175+176+177+178+179+180+181+182+183+184+185+201+203+204+205+206+207+208+209+210+211+212
color ucsf_pink1, in2to4int
### in 5 or more interfaces
select morethan4, 3gj0 and resi 12+19+29+35+37+38+39+40+41+43+44+45+46+47+62+64+68+69+70+71+72+74+75+76+77+78+79+81+82+91+93+94+95+96+99+100+103+106+107+110+111+113+126+127+128+129+130+134+137+139+140+141+142+143+145+146+147+152+153+154+155+156+159+163+166+173
color ucsf_purple1, morethan4
hide all
show surface
set_view (\
    -0.665958822,    0.680874586,   -0.304802597,\
    -0.669537663,   -0.365375221,    0.646693468,\
     0.328949124,    0.634752274,    0.699200273,\
     0.000000000,    0.000000000, -163.817398071,\
    17.855573654,    1.621843338,  -20.121271133,\
   129.154876709,  198.479919434,  -20.000000000 )
show sticks, not polymer
color black, not polymer
set valence, off
set depth_cue, off
set specular, off
set ray_shadows, off
set surface_quality, 2
png /Users/tperica/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_interface_patches/all_interface_residue_freq.png, dpi = 1200, 0, 0, -1
rotate x, 90
png /Users/tperica/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_interface_patches/all_interface_residue_freq_90deg.png, dpi = 1200, 0, 0, -1
rotate x, 90
png /Users/tperica/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_interface_patches/all_interface_residue_freq_180deg.png, dpi = 1200, 0, 0, -1
rotate x, 90
png /Users/tperica/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_interface_patches/all_interface_residue_freq_270deg.png, dpi = 1200, 0, 0, -1


