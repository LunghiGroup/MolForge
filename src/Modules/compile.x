#! /bin/bash

 mpif90 -g -O3 -c def_lists.f90 
 mpif90 -g -O3 -c def_random_numbers.f90
 mpif90 -g -O3 -c def_angmom.f90  
 mpif90 -g -O3 -c def_general_types.f90  
 mpif90 -g -O3 -c def_mpi.f90      
 mpif90 -g -O3 -c def_blacs.f90   
 mpif90 -g -O3 -c def_lapack.f90         
 mpif90 -g -O3 -c def_orthonorm.f90         
 mpif90 -g -O3 -c def_scalapack.f90  
 mpif90 -g -O3 -c def_utils.f90
 mpif90 -g -O3 -c def_dist.f90
 mpif90 -g -O3 -c def_parser.f90   
 mpif90 -g -O3 -c def_sparse.f90      
 mpif90 -g -O3 -c def_stevens.f90
 mpif90 -g -O3 -c def_rotations.f90  
 mpif90 -g -O3 -c def_nets.f90
 mpif90 -g -O3 -c def_bispectrum.f90
 mpif90 -g -O3 -c def_descriptors.f90
 mpif90 -g -O3 -c def_forcefields.f90
 mpif90 -g -O3 -c def_lattice.f90        
 mpif90 -g -O3 -c def_atoms.f90   
 mpif90 -g -O3 -c def_proj_disp.f90      

 mpif90 -g -O3 -c def_exponential.f90
 mpif90 -g -O3 -c def_spinham.f90     
 mpif90 -g -O3 -c def_phonons.f90  
 mpif90 -g -O3 -c def_spinphonon.f90  
 mpif90 -g -O3 -c def_spins_dist_rs.f90  
 mpif90 -fmax-errors=10 -g -O3 -c def_hilbert_dist.f90  

 mpif90 -g -O3 -c def_function.f90
 mpif90 -g -O3 -c def_grad.f90
 mpif90 -g -O3 -c def_swarm.f90
 mpif90 -g -O3 -c def_ff_training.f90

 mpif90 -g -O3 -c def_spinham_map.f90  
 mpif90 -g -O3 -c def_pulses.f90     

 ar crv libmolforge.a *.o

 exit

 mpif90 -g -O3 -c def_forcefields.f90
 mpif90 -g -O3 -c def_atoms.f90   
 mpif90 -g -O3 -c def_ff_training.f90

 ar crv libmolforge.a *.o

 exit


 for name in def_lists.f90 def_random_numbers.f90 def_angmom.f90 def_general_types.f90 def_mpi.f90 def_blacs.f90 def_lapack.f90 def_orthonorm.f90 def_scalapack.f90 def_utils.f90 def_dist.f90 def_parser.f90 def_sparse.f90 def_stevens.f90 def_rotations.f90 def_nets.f90 def_bispectrum.f90 def_descriptors.f90 def_forcefields.f90 def_lattice.f90 def_atoms.f90 def_proj_disp.f90 def_exponential.f90 def_spinham.f90 def_phonons.f90 def_spinphonon.f90 def_spins_dist_rs.f90 def_hilbert_dist.f90 def_function.f90 def_grad.f90 def_swarm.f90 def_ff_training.f90 def_spinham_map.f90 def_pulses.f90
 do

##  export uptodate=$(diff ${name} last_update/${name} | wc -l | awk '{print $1}')
##  if [ ${uptodate} -ne 0 ]   
##  then
   echo "mpif90 -fmax-errors=10 -g -O3 -c "${name}
   mpif90 -fmax-errors=10 -g -O3 -c ${name} 
   export result=$?
   if [ ${result} -eq 0 ]   
   then
    echo "cp "${name}" last_updata/"
    cp ${name} last_update/
   fi
##  fi

 done
  
 ar crv libmolforge.a *.o

 exit

