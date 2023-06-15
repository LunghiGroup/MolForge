.SUFFIXES:

SOURCES = def_lammps.f90 main.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

Snap : ${OBJECTS}
	$(FC) $(SNAP_DIR)/$(OBJECTS) $(SNAP_INC_FLAG) $(INC_FLAG) $(SNAP_LINK_FLAG) $(LINK_FLAG) -o Snap.x
	cp $(SNAP_DIR)/Snap.x $(EXE_DIR)/

def_lammps.o : def_lammps.f90
	$(FC) $(FFLAGS) -c $(SNAP_DIR)/$(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(SNAP_INC_FLAG)

main.o : def_lammps.o main.f90
	$(FC) $(FFLAGS) -c $(SNAP_DIR)/$(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(SNAP_INC_FLAG)

def_lammps.f90:

main.f90: def_lammps.o
