.SUFFIXES:

export FC = mpif90
export FFLAGS = -g -O0 -fcheck=all -fbacktrace

export LIBS = -lscalapack-openmpi -llapack
export INC_FLAG =## -I$(INC_DIR) 
export LINK_FLAG = $(LIBS) ## -L${LIB_DIR} ${LIBS}

## Set this keyword to "MKL" if you are using MKL scalapack instead of GNU libs.
export SCALAPACK_FLAG =

######## The rest of the flags should not be touched

export MOLFORGE_DIR = $(shell pwd )

export MODULES_DIR = $(MOLFORGE_DIR)/Modules

export SPIRAL_DIR = $(MOLFORGE_DIR)/Spiral
export PHONDY_DIR = $(MOLFORGE_DIR)/PhonDy
export INTENSOR_DIR = $(MOLFORGE_DIR)/Intensor
export TOOLS_DIR = $(MOLFORGE_DIR)/Tools
export SNAP_DIR =$(MOLFORGE_DIR)/Snap

export LAMMPS_DIR = /home/valerio/Desktop/lammps
export LAMMPS_F90_DIR = $(LAMMPS_DIR)/examples/COUPLE/fortran2
export LAMMPS_MPI_DIR = $(LAMMPS_DIR)/src

export EXE_DIR = $(MOLFORGE_DIR)/bin

export MOLFORGE_INC_DIR = $(MODULES_DIR)
export MOLFORGE_LIB_DIR = $(MOLFORGE_DIR)/libs
export MOLFORGE_LIBS = -lMolForge 

export MOLFORGE_INC_FLAG = -I$(MOLFORGE_INC_DIR) 
export MOLFORGE_LINK_FLAG = -L$(MOLFORGE_LIB_DIR) $(MOLFORGE_LIBS)

export MPI_LIB = /usr/lib64/openmpi/lib
export SNAP_LIBS = -lblas -llammps_fortran -llammps_mpi -lMolForge -llapack -lmpi_cxx -lstdc++ -lm -ldftd3
export DFTD3_LIB= /home/valerio/Desktop/dftd3-lib-0.9/lib

export SNAP_INC_FLAG = -I$(MODULES_DIR) -I$(LAMMPS_F90_DIR) -I$(DFTD3_LIB)
export SNAP_LINK_FLAG = -L$(MPI_LIB) -L$(LAMMPS_MPI_DIR) -L$(LAMMPS_F90_DIR) -L$(MOLFORGE_LIB_DIR) -L$(DFTD3_LIB) $(SNAP_LIBS)




all : libMolForge Spiral PhonDy Tools

libMolForge :
	mkdir -p $(MOLFORGE_DIR)/libs
	cd $(MODULES_DIR) && $(MAKE) $@

Spiral : libMolForge
	mkdir -p $(MOLFORGE_DIR)/bin
	cd $(SPIRAL_DIR) && $(MAKE) $@

PhonDy : libMolForge
	mkdir -p $(MOLFORGE_DIR)/bin
	cd $(PHONDY_DIR) && $(MAKE) $@

Tools : libMolForge
	mkdir -p $(MOLFORGE_DIR)/bin
	cd $(TOOLS_DIR) && $(MAKE) all

Snap : libMolForge
	mkdir -p $(MOLFORGE_DIR)/bin
	cd $(SNAP_DIR) && $(MAKE) $@

clean : clean_MolForge clean_Spiral clean_PhonDy clean_Tools 

clean_MolForge :
	rm -f $(MODULES_DIR)/*.o
	rm -f $(MODULES_DIR)/*.mod
	rm -f $(MOLFORGE_LIB_DIR)/libMolForge.a

clean_Spiral : 
	rm -f $(SPIRAL_DIR)/*.o
	rm -f $(SPIRAL_DIR)/*.mod
	rm -f $(EXE_DIR)/Spiral.x

clean_PhonDy : 
	rm -f $(PHONDY_DIR)/*.o
	rm -f $(PHONDY_DIR)/*.mod
	rm -f $(EXE_DIR)/PhonDy.x

clean_Snap :
	rm -f $(SNAP_DIR)/*.o
	rm -f $(SNAP_DIR)/*.mod
	rm -f $(EXE_DIR)/Snap.x

clean_Tools :
	cd $(TOOLS_DIR) && $(MAKE) clean
