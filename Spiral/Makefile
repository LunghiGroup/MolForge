.SUFFIXES:

SOURCES = variables.f90 parse_input.f90 Spiral.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

Spiral : ${OBJECTS}
	$(shell cp Spiral.x ../bin/)
	$(FC) $(SPIRAL_DIR)/$(OBJECTS) -I$(MOLFORGE_INC_DIR) -L$(MOLFORGE_LIB_DIR) $(MOLFORGE_LIBS) $(FFLAGS) -o Spiral.x


variables.o : variables.f90
	$(FC) -I$(MOLFORGE_INC_DIR) $(FFLAGS) -c $(SPIRAL_DIR)/$(subst .o,.f90,$@) 

variables.f90 : 


parse_input.o : variables.o parse_input.f90
	$(FC) -I$(MOLFORGE_INC_DIR) $(FFLAGS) -c $(SPIRAL_DIR)/$(subst .o,.f90,$@) 

parse_input.f90 : variables.o


Spiral.o : variables.o parse_input.o Spiral.f90
	$(FC) -I$(MOLFORGE_INC_DIR)  $(FFLAGS) -c $(SPIRAL_DIR)/$(subst .o,.f90,$@) 

Spiral.f90 : variables.o parse_input.o



