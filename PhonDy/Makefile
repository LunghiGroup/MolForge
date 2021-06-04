.SUFFIXES:

SOURCES = PhonDy.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

PhonDy : ${OBJECTS}
	$(FC) $(PHONDY_DIR)/$(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_DIR) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o PhonDy.x
	mv $(PHONDY_DIR)/PhonDy.x $(EXE_DIR)/

PhonDy.o : PhonDy.f90
	$(FC) $(FFLAGS) -c $(PHONDY_DIR)/$(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG) 

PhonDy.f90 : 



