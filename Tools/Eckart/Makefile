.SUFFIXES:

SOURCES = Eckart.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

Eckart : $(OBJECTS)
	$(FC) $(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_FLAG) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o Eckart.x
	mv Eckart.x $(EXE_DIR)/

Eckart.o : Eckart.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

Eckart.f90 : 

clean :
	rm -f *.o
	rm -f *.mod
	rm -f $(EXE_DIR)/Eckart.x
