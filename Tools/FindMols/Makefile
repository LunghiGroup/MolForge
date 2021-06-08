.SUFFIXES:

SOURCES = find_mols.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

FindMols : $(OBJECTS)
	$(FC) $(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_FLAG) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o find_mols.x
	mv find_mols.x $(EXE_DIR)/

find_mols.o : find_mols.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

find_mols.f90 : 

clean :
	rm -f *.o
	rm -f *.mod
	rm -f $(EXE_DIR)/find_mols.x
