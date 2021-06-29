.SUFFIXES:

SOURCES = getFC.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

getFC : $(OBJECTS)
	$(FC) $(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_FLAG) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o getFC.x
	mv getFC.x $(EXE_DIR)/

getFC.o : getFC.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

getFC.f90 : 

clean :
	rm -f *.o
	rm -f *.mod
	rm -f $(EXE_DIR)/getFC.x
