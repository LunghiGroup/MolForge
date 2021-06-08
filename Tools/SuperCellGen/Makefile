.SUFFIXES:

SOURCES = pbc.f90 SuperCellGen.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

SuperCellGen : $(OBJECTS)
	$(FC) $(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_FLAG) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o SuperCellGen.x
	mv SuperCellGen.x $(EXE_DIR)/

SuperCellGen.o : pbc.o SuperCellGen.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

SuperCellGen.f90 : 

pbc.o : pbc.f90 
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

pbc.f90 : 

clean :
	rm -f *.o
	rm -f *.mod
	rm -f $(EXE_DIR)/SuperCellGen.x
