.SUFFIXES:

ifeq ($(SCALAPACK_FLAG),MKL)
$(info $$SCALAPACK_FLAG is [${SCALAPACK_FLAG}])
SOURCES = def_lists.f90 def_random_numbers.f90 def_angmom.f90 def_general_types.f90 def_mpi.f90 def_blacs.f90 def_lapack.f90 def_orthonorm.f90  def_scalapack_mkl.f90 def_utils.f90 def_dist.f90 def_parser.f90 def_sparse.f90 def_stevens.f90 def_rotations.f90 def_nets.f90 def_bispectrum.f90 def_descriptors.f90 def_forcefields.f90 def_lattice.f90 def_atoms.f90 def_proj_disp.f90 def_exponential.f90 def_spinham.f90 def_phonons.f90 def_spinphonon.f90 def_spins_dist_rs.f90 def_hilbert_dist.f90 def_function.f90 def_grad.f90 def_swarm.f90 def_ff_training.f90 def_spinham_map.f90 def_pulses.f90 def_mlmodel.f90 def_mlmodel_training.f90 def_sh_mlmodel_training.f90
else
$(info $$SCALAPACK_FLAG is [${SCALAPACK_FLAG}]) 
SOURCES = def_lists.f90 def_random_numbers.f90 def_angmom.f90 def_general_types.f90 def_mpi.f90 def_blacs.f90 def_lapack.f90 def_orthonorm.f90  def_scalapack.f90 def_utils.f90 def_dist.f90 def_parser.f90 def_sparse.f90 def_stevens.f90 def_rotations.f90 def_nets.f90 def_bispectrum.f90 def_descriptors.f90 def_forcefields.f90 def_lattice.f90 def_atoms.f90 def_proj_disp.f90 def_exponential.f90 def_spinham.f90 def_phonons.f90 def_spinphonon.f90 def_spins_dist_rs.f90 def_hilbert_dist.f90 def_function.f90 def_grad.f90 def_swarm.f90 def_ff_training.f90 def_spinham_map.f90 def_pulses.f90 def_mlmodel.f90 def_mlmodel_training.f90 def_sh_mlmodel_training.f90
endif

OBJECTS = $(subst .f90,.o,$(SOURCES))

libMolForge : $(OBJECTS)
	ar crv $@.a $(OBJECTS)
	mv $(MODULES_DIR)/libMolForge.a $(MOLFORGE_LIB_DIR)/

########## Make object files

def_lists.o : def_lists.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_random_numbers.o : def_random_numbers.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_angmom.o : def_angmom.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)
	
def_general_types.o : def_general_types.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_mpi.o : def_mpi.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_lapack.o : def_lapack.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_blacs.o : def_mpi.o def_blacs.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_orthonorm.o : def_lapack.o def_orthonorm.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

ifeq ($(SCALAPACK_FLAG),MKL)
def_scalapack_mkl.o : def_mpi.o def_blacs.o def_scalapack_mkl.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)
else
def_scalapack.o : def_mpi.o def_blacs.o def_scalapack.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)
endif

def_utils.o : def_utils.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_dist.o : def_utils.o def_dist.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_parser.o : def_parser.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_sparse.o : def_lists.o def_mpi.o def_blacs.o def_sparse.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_stevens.o : def_angmom.o def_stevens.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_rotations.o: def_angmom.o def_rotations.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_nets.o : def_nets.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_bispectrum.o : def_angmom.o def_bispectrum.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_descriptors.o : def_bispectrum.o def_lapack.o def_bispectrum.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_forcefields.o : def_nets.o def_descriptors.o def_forcefields.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_mlmodel.o : def_nets.o def_descriptors.o def_mlmodel.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_lattice.o : def_mpi.o def_lapack.o def_lattice.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_atoms.o : def_lattice.o def_descriptors.o def_forcefields.o def_lists.o def_sparse.o def_mpi.o def_utils.o def_general_types.o def_atoms.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_lammps.o : def_atoms.o def_lammps.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_proj_disp.o : def_lapack.o def_proj_disp.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_exponential.o : def_lists.o def_mpi.o def_blacs.o def_stevens.o def_rotations.o def_lapack.o def_utils.o def_exponential.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_spinham.o : def_lists.o def_mpi.o def_blacs.o def_utils.o def_stevens.o def_rotations.o def_lapack.o def_spinham.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_phonons.o : def_mpi.o def_lapack.o def_utils.o def_atoms.o def_lattice.o def_proj_disp.o def_phonons.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_spinphonon.o : def_lists.o def_mpi.o def_dist.o def_blacs.o def_spinham.o def_stevens.o def_phonons.o def_proj_disp.o def_spinphonon.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_spins_dist_rs.o : def_mpi.o def_blacs.o def_lattice.o def_spinham.o def_spinphonon.o def_spins_dist_rs.f90 
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

ifeq ($(SCALAPACK_FLAG),MKL)
def_hilbert_dist.o : def_mpi.o def_blacs.o def_lapack.o def_scalapack_mkl.o def_lists.o def_sparse.o def_dist.o def_utils.o def_lattice.o def_atoms.o def_phonons.o def_spinham.o def_spinphonon.o def_spins_dist_rs.o def_exponential.o def_pulses.o def_hilbert_dist.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)
else
def_hilbert_dist.o : def_mpi.o def_blacs.o def_lapack.o def_scalapack.o def_lists.o def_sparse.o def_dist.o def_utils.o def_lattice.o def_atoms.o def_phonons.o def_spinham.o def_spinphonon.o def_spins_dist_rs.o def_exponential.o def_pulses.o def_hilbert_dist.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)
endif

def_function.o : def_function.f90 
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_grad.o : def_function.o def_random_numbers.o def_grad.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_swarm.o : def_function.o def_random_numbers.o def_swarm.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_ff_training.o : def_nets.o def_function.o def_atoms.o def_general_types.o def_random_numbers.o def_forcefields.o def_ff_training.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_mlmodel_training.o : def_nets.o def_function.o def_atoms.o def_general_types.o def_random_numbers.o def_mlmodel.o def_descriptors.o def_mlmodel_training.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_sh_mlmodel_training.o : def_nets.o def_function.o def_atoms.o def_general_types.o def_random_numbers.o def_descriptors.o def_mlmodel.o def_mlmodel_training.o def_sh_mlmodel_training.f90 
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_spinham_map.o : def_lapack.o def_angmom.o def_stevens.o def_swarm.o def_orthonorm.o def_hilbert_dist.o def_spinham_map.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)

def_pulses.o : def_mpi.o def_blacs.o def_lists.o def_pulses.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(INC_FLAG)


########## Make source fiels


def_lists.f90 : 

def_random_numbers.f90 : 

def_angmom.f90 : 
	
def_general_types.f90 : 

def_mpi.f90 : 

def_lapack.f90 : 

def_blacs.f90 : def_mpi.o 

def_orthonorm.f90 : def_lapack.o 

def_scalapack.f90 : def_mpi.o def_blacs.o 

def_scalapack_mkl.f90 : def_mpi.o def_blacs.o 

def_utils.f90 : 

def_dist.f90 : def_utils.o 

def_parser.f90 : 

def_sparse.f90 : def_lists.o def_mpi.o def_blacs.o

def_stevens.f90 : def_angmom.o 

def_rotations.f90 : def_angmom.o

def_nets.f90 :

def_bispectrum.f90 : def_angmom.o 

def_descriptors.f90 : def_bispectrum.o def_lapack.o 

def_forcefields.f90 : def_nets.o def_descriptors.o 

def_mlmodel.f90 : def_nets.o def_descriptors.o 

def_lattice.f90 : def_mpi.o def_lapack.o 

def_atoms.f90 : def_lattice.o def_descriptors.o def_forcefields.o def_lists.o def_sparse.o def_mpi.o def_utils.o def_general_types.o  def_mlmodel.o

def_proj_disp.f90 : def_lapack.o 

def_exponential.f90 : def_lists.o def_mpi.o def_blacs.o def_stevens.o def_rotations.o def_lapack.o def_utils.o 

def_spinham.f90 : def_lists.o def_mpi.o def_blacs.o def_utils.o def_stevens.o def_rotations.o def_lapack.o 

def_phonons.f90 : def_mpi.o def_lapack.o def_utils.o def_atoms.o def_lattice.o def_proj_disp.o 

def_spinphonon.f90 : def_lists.o def_mpi.o def_dist.o def_blacs.o def_spinham.o def_stevens.o def_phonons.o def_proj_disp.o 

def_spins_dist_rs.f90 : def_mpi.o def_blacs.o def_lattice.o def_spinham.o def_spinphonon.o 

ifeq ($(SCALAPACK_FLAG),MKL)
def_hilbert_dist.f90 : def_mpi.o def_blacs.o def_lapack.o def_scalapack_mkl.o def_lists.o def_sparse.o def_dist.o def_utils.o def_lattice.o def_atoms.o def_phonons.o def_spinham.o def_spinphonon.o def_spins_dist_rs.o def_exponential.o def_pulses.o 
else
def_hilbert_dist.f90 : def_mpi.o def_blacs.o def_lapack.o def_scalapack.o def_lists.o def_sparse.o def_dist.o def_utils.o def_lattice.o def_atoms.o def_phonons.o def_spinham.o def_spinphonon.o def_spins_dist_rs.o def_exponential.o def_pulses.o 
endif

def_function.f90 : 

def_grad.f90 : def_function.o def_random_numbers.o

def_swarm.f90 : def_function.o def_random_numbers.o 

def_ff_training.f90 : def_nets.o def_function.o def_atoms.o def_general_types.o def_random_numbers.o def_forcefields.o

def_mlmodel_training.f90 : def_nets.o def_function.o def_atoms.o def_general_types.o def_descriptors.o def_random_numbers.o def_mlmodel.o

def_sh_mlmodel_training.f90 : def_nets.o def_function.o def_atoms.o def_general_types.o def_random_numbers.o def_descriptors.o def_mlmodel.o def_mlmodel_training.o 

def_spinham_map.f90 : def_lapack.o def_angmom.o def_stevens.o def_swarm.o def_orthonorm.o def_hilbert_dist.o 

def_pulses.f90 : def_mpi.o def_blacs.o def_lists.o 

