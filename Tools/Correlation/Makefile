.SUFFIXES:

SOURCES = Correlation.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

Correlation : $(OBJECTS)
	$(FC) $(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_FLAG) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o Correlation.x
	mv Correlation.x $(EXE_DIR)/

Correlation.o : Correlation.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

Correlation.f90 : 

clean :
	rm -f *.o
	rm -f *.mod
	rm -f $(EXE_DIR)/Correlation.x
