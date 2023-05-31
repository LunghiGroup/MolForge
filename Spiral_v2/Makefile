.SUFFIXES:

SOURCES = Hdyn.f90 new_hilbert.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

Hdyn : ${OBJECTS}
	$(FC) $(HDYN_DIR)/$(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_FLAG) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o Hdyn.x
	mv $(HDYN_DIR)/Hdyn.x $(EXE_DIR)/

new_hilbert.o : new_hilbert.f90
	$(FC) $(FFLAGS) -c $(SPIRAL_DIR)/$(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG) 

new_hilbert.f90 : 

Hdyn.o : new_hilbert.o Hdyn.f90
	$(FC) $(FFLAGS) -c $(SPIRAL_DIR)/$(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG) 

Hdyn.f90 : new_hilbert.o



