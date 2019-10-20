# get two hemoglobin beta chains, the "ligand" will be the heme hetero group
fetch 1hbb, async=0
create conf1, chain B
create conf2, chain D
delete 1hbb

# important: identifiers must be the same
alter conf2, chain="B"
alter all, segi=""

# optional: styling
as cartoon
show sticks, not polymer
show nb_spheres

# superpose and morph
align conf1, conf2
#morph mout, conf1, conf2, match=in, refinement=0