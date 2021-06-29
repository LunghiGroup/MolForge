.SUFFIXES:

SOURCES = fitcharges.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

fitcharges : $(OBJECTS)
	$(FC) $(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_FLAG) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o fitcharges.x
	mv fitcharges.x $(EXE_DIR)/

fitcharges.o : fitcharges.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

fitcharges.f90 : 

clean :
	rm -f *.o
	rm -f *.mod
	rm -f $(EXE_DIR)/fitcharges.x
