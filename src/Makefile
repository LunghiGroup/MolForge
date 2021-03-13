
.SUFFIXES:

MOLFORGE_DIR = /home/Alessandro/ERC/MolForge
SRC_DIR = $(MOLFORGE_DIR)/src
MODULES_DIR = $(SRC_DIR)/Modules
OBJECTS_DIR = $(SRC_DIR)/Objects
INTERFACES_DIR = $(SRC_DIR)/Interfaces

FC = mpif90
FC_FLAGS = -g -fmax-errors=10 -O3

SCALAPACK_FLAG = GNU ## Has to be "MKL" or "GNU", any other value will default to GNU

SOURCES = def_lists.f90 def_random_numbers.f90 # def_angmom.f90 def_general_types.f90 def_mpi.f90 def_blacs.f90 def_lapack.f90 def_orthonorm.f90 def_scalapack.f90 def_utils.f90 def_dist.f90 def_parser.f90 def_sparse.f90 def_stevens.f90 def_rotations.f90 def_nets.f90 def_bispectrum.f90 def_descriptors.f90 def_forcefields.f90 def_lattice.f90 def_atoms.f90 def_proj_disp.f90 def_exponential.f90 def_spinham.f90 def_phonons.f90 def_spinphonon.f90 def_spins_dist_rs.f90 def_hilbert_dist.f90 def_function.f90 def_grad.f90 def_swarm.f90 def_ff_training.f90 def_spinham_map.f90 def_pulses.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

libMolForge.a : $(OBJECTS)
	ar crv $@ $^
##	@mkdir -p $(OBJECTS_DIR)
##	@mkdir -p $(INTERFACES_DIR)
##	@mv $^ $(OBJECTS_DIR)/
##	@mv *.mod $(INTERFACES_DIR)/

def_lists.o : 
	$(FC) $(FC_FLAGS) -c $(MODULES_DIR)/$(subst .o,.f90,$@)

def_random_numbers.o : 
	$(FC) $(FC_FLAGS) -c $(MODULES_DIR)/$(subst .o,.f90,$@)

