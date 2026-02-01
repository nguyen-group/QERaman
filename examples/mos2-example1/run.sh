#mpirun -np 48 pw.x <scf.in >scf.out
#mpirun -np 48 bands_mat.x <bands_mat.in> bands_mat.out
mpirun -np 48 ph_mat.x <ph_mat.in> ph_mat.out
raman.x <raman.in> raman.out
