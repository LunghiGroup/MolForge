#! /bin/bash

 mpif90 -g -O2 -c def_lists.f90 
 mpif90 -g -O2 -c def_random_numbers.f90
 mpif90 -g -O2 -c def_angmom.f90  
 mpif90 -g -O2 -c def_general_types.f90  
 mpif90 -g -O2 -c def_mpi.f90      
 mpif90 -g -O2 -c def_blacs.f90   
 mpif90 -g -O2 -c def_lapack.f90         
 mpif90 -g -O2 -c def_orthonorm.f90         
 mpif90 -g -O2 -c def_scalapack.f90  
 mpif90 -g -O2 -c def_utils.f90
 mpif90 -g -O2 -c def_dist.f90
 mpif90 -g -O2 -c def_parser.f90   
 mpif90 -g -O2 -c def_sparse.f90      
 mpif90 -g -O2 -c def_stevens.f90
 mpif90 -g -O2 -c def_rotations.f90  
 mpif90 -g -O2 -c def_nets.f90
 mpif90 -g -O2 -c def_bispectrum.f90
 mpif90 -g -O2 -c def_descriptors.f90
 mpif90 -g -O2 -c def_forcefields.f90
 mpif90 -g -O2 -c def_lattice.f90        
 mpif90 -g -O2 -c def_atoms.f90   
 mpif90 -g -O2 -fcheck=all -c def_proj_disp.f90      

 mpif90 -g -O2 -c def_pulses.f90     
 mpif90 -g -O2 -c def_exponential.f90
 mpif90 -g -O2 -c def_spinham.f90     
 mpif90 -g -O2 -c def_phonons.f90  
 mpif90 -g -O2 -c def_spinphonon.f90  
 mpif90 -g -O2 -c def_spins_dist_rs.f90  
 mpif90 -fmax-errors=10 -g -O2 -c def_hilbert_dist.f90  

 mpif90 -g -O2 -c def_function.f90
 mpif90 -g -O2 -c def_grad.f90
 mpif90 -g -O2 -c def_swarm.f90
 mpif90 -g -O2 -c def_ff_training.f90
 mpif90 -g -O2 -c def_spinham_map.f90  

 ar crv libmolforge.a *.o

