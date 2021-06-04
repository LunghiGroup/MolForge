.SUFFIXES:

SOURCES = variables.f90 parse_input.f90 Spiral.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

Spiral : ${OBJECTS}
	$(FC) $(SPIRAL_DIR)/$(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_FLAG) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o Spiral.x
	mv $(SPIRAL_DIR)/Spiral.x $(EXE_DIR)/


variables.o : variables.f90
	$(FC) $(FFLAGS) -c $(SPIRAL_DIR)/$(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

variables.f90 : 


parse_input.o : variables.o parse_input.f90
	$(FC) $(FFLAGS) -c $(SPIRAL_DIR)/$(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG) 

parse_input.f90 : variables.o


Spiral.o : variables.o parse_input.o Spiral.f90
	$(FC) $(FFLAGS) -c $(SPIRAL_DIR)/$(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG) 

Spiral.f90 : variables.o parse_input.o



