# Making microtubule

# 14 PF

# Fetch 4hna to 14PF microtubule map (emd-8998), fit 4hna in PF1 next to the seam, save the file as another pdb (4hna_align_PF1.pdb). Close and load the 4hna_align_PF1.pdb again.
# If the map is 512^3 with pixel size 1 => center in Z = 256

# Symmetrization
sym #1 group h,8.72,-25.77,14*shift,3,82.0 surf true axis z center 256,256,256 res 5

# Making a map of kinesin only
select #2
# Action/Ribbon Show
molmap :*.K 15

# Resample on the map
vop resample #3 onGrid #0

# Save as map

# 13 PF

