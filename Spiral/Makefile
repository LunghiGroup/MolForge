.SUFFIXES:

SOURCES = variables.f90 parse_input.f90 Spiral.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

Spiral : ${OBJECTS}
	$(FC) $(LIBS) $(FFLAGS) $(SPIRAL_DIR)/$(OBJECTS) -o Spiral.x


variables.o : variables.f90
	$(FC) $(FFLAGS) -c $(MODULES_DIR)/$(subst .o,.f90,$@) 

variables.f90 : 


parse_input.o : variables.f90 parse_input.f90
	$(FC) $(FFLAGS) -c $(MODULES_DIR)/$(subst .o,.f90,$@) 

parse_input.f90 : 


Spiral.o : variables.f90 parse_input.f90 Spiral.f90
	$(FC) $(FFLAGS) -c $(MODULES_DIR)/$(subst .o,.f90,$@) 

Spiral.f90 : 



