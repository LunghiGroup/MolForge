.SUFFIXES:

SOURCES = CF.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

CF : $(OBJECTS)
	$(FC) $(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_FLAG) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o CF.x
	mv CF.x $(EXE_DIR)/

CF.o : CF.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

CF.f90 : 

clean :
	rm -f *.o
	rm -f *.mod
	rm -f $(EXE_DIR)/CF.x
