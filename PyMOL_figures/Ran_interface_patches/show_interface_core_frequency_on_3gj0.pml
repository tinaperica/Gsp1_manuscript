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
select zerointerfaces, 3gj0 and resi 1+2+5+6+7+8+9+10+11+12+13+14+15+16+17+18+21+22+23+24+25+26+27+28+29+30+36+37+38+40+42+48+49+50+52+54+58+59+60+61+62+63+65+66+68+79+80+82+83+84+85+86+87+88+89+90+91+92+97+98+101+104+105+108+109+112+113+114+115+116+117+118+119+120+121+122+124+127+131+134+136+138+144+146+147+148+149+150+151+156+160+161+162+163+164+165+167+170+172+173+174+175+182+183+185+186+187+188+189+190+191+192+193+194+195+196+197+198+199+200+202+203+206+207+208+213+214+215+216+217+218+219
color white, zerointerfaces
### in 1 interface
select oneinterface, 3gj0 and resi 31+32+41+44+45+46+51+53+100+128+132+142+143+157+171+177+178+201+210+211+212
color ucsf_pink3, oneinterface
### in 2-4 interfaces
select in2to4int, 3gj0 and resi 107+154+158+166+169+20+33+64+72+73+95+96+102+103+126+130+133+137+141+155+3+4+34+39+43+55+56+57+67+69+70+71+93+99+123+125+129+135+139+152+153+168+176+179+180+181+184+204+205+209
color ucsf_pink1, in2to4int
### in 5 or more interfaces
select morethan4, 3gj0 and resi 78+74+75+81+110+47+111+77+19+35+76+94+106+140+145+159
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
png /Users/tperica/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_interface_patches/core_residue_freq.png, dpi = 1200, 0, 0, -1
rotate x, 90
png /Users/tperica/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_interface_patches/core_residue_freq_90deg.png, dpi = 1200, 0, 0, -1
rotate x, 90
png /Users/tperica/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_interface_patches/core_residue_freq_180deg.png, dpi = 1200, 0, 0, -1
rotate x, 90
png /Users/tperica/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_interface_patches/core_residue_freq_270deg.png, dpi = 1200, 0, 0, -1



png /Users/tperica/Box Sync/kortemmelab/home/tina/Gsp1_manuscript/PyMOL_figures/Ran_interface_patches/cartoon.png, dpi = 1200, 0, 0, -1
