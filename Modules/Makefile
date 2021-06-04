.SUFFIXES:

SOURCES = def_lists.f90 def_random_numbers.f90 def_angmom.f90 def_general_types.f90 def_mpi.f90 def_blacs.f90 def_lapack.f90 def_orthonorm.f90 ## def_scalapack.f90 def_utils.f90 def_dist.f90 def_parser.f90 def_sparse.f90 def_stevens.f90 def_rotations.f90 def_nets.f90 def_bispectrum.f90 def_descriptors.f90 def_forcefields.f90 def_lattice.f90 def_atoms.f90 def_proj_disp.f90 def_exponential.f90 def_spinham.f90 def_phonons.f90 def_spinphonon.f90 def_spins_dist_rs.f90 def_hilbert_dist.f90 def_function.f90 def_grad.f90 def_swarm.f90 def_ff_training.f90 def_spinham_map.f90 def_pulses.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

libMolForge : $(OBJECTS)
	ar crv $@.a $(OBJECTS)


########## Make object files


def_lists.o : def_lists.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) 

def_random_numbers.o : def_random_numbers.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@)

def_angmom.o : def_angmom.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@)
	
def_general_types.o : def_general_types.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@)

def_mpi.o : def_mpi.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@)

def_blacs.o : def_blacs.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@)

def_lapack.o : def_lapack.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@)

def_orthonorm.o : def_lapack.o def_orthonorm.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@)


########## Make source fiels


def_lists.f90 : 

def_random_numbers.f90 : 

def_angmom.f90 :
	
def_general_types.f90 : 

def_mpi.f90 :

def_blacs.f90 : 

def_lapack.f90 :

def_orthonorm.f90 : def_lapack.o
