.SUFFIXES:

SOURCES = Intensor.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

Intensor : ${OBJECTS}
	$(FC) $(INTENSOR_DIR)/$(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_FLAG) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o Intensor.x
	mv $(INTENSOR_DIR)/Intensor.x $(EXE_DIR)/


Intensor.o : Intensor.f90
	$(FC) $(FFLAGS) -c $(INTENSOR_DIR)/$(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

Intensor.f90 : 


