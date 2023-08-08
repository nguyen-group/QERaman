rm -rf ../tmp
mpirun -np 48 pw.x <scf.in >scf.out
mpirun -np 48 bands_mat.x <bands.in> bands.out
mpirun -np 48 ph_mat.x <ph.in> ph.out
raman.x <raman.in> raman.out
