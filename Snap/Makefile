.SUFFIXES:

SOURCES = def_kind.f90 def_parameters.f90 def_lammps.f90 def_potential.f90 def_max.f90 def_VdW.f90 def_SNAP_fit.f90 def_SNAP_FF.f90 def_molecular_dynamics.f90 main.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

Snap : ${OBJECTS}
	$(FC) -g $(SNAP_DIR)/$(OBJECTS) $(SNAP_INC_FLAG) $(INC_FLAG) $(SNAP_LINK_FLAG) $(LINK_FLAG) -o Snap.x
	cp $(SNAP_DIR)/Snap.x $(EXE_DIR)/

def_kind.o : def_kind.f90
	$(FC) $(FFLAGS) -c $(SNAP_DIR)/$(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(SNAP_INC_FLAG)

def_parameters.o: def_parameters.f90
	$(FC) $(FFLAGS) -c $(SNAP_DIR)/$(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(SNAP_INC_FLAG)

def_lammps.o : def_lammps.f90
	$(FC) $(FFLAGS) -c $(SNAP_DIR)/$(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(SNAP_INC_FLAG)

def_potential.o: def_potential.f90
	$(FC) $(FFLAGS) -c $(SNAP_DIR)/$(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(SNAP_INC_FLAG)

def_max.o : def_max.f90
	$(FC) $(FFLAGS) -c $(SNAP_DIR)/$(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(SNAP_INC_FLAG)

def_SNAP_fit.o : def_SNAP_fit.f90
	$(FC) $(FFLAGS) -c $(SNAP_DIR)/$(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(SNAP_INC_FLAG)

def_SNAP_FF.o : def_SNAP_FF.f90
	$(FC) $(FFLAGS) -c $(SNAP_DIR)/$(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(SNAP_INC_FLAG)

def_VdW.o: def_VdW.f90
	$(FC) $(FFLAGS) -c $(SNAP_DIR)/$(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(SNAP_INC_FLAG)

def_molecular_dynamics.o: def_molecular_dynamics.f90
	$(FC) $(FFLAGS) -c $(SNAP_DIR)/$(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(SNAP_INC_FLAG)

main.o : main.f90
	$(FC) $(FFLAGS) -c $(SNAP_DIR)/$(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(SNAP_INC_FLAG)

def_rescale.f90:

def_kind.f90:

def_lammps.f90:

def_max.f90:

def_parameters.f90:

def_SNAP_fit.f90:

def_VdW.f90:

def_SNAP_FF.f90:
	
def_potential.f90:

def_molecular_dynamics.f90:

main.f90: def_kind.o def_lammps.o

