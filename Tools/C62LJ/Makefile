.SUFFIXES:

SOURCES = C62LJ.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

C62LJ : $(OBJECTS)
	$(FC) $(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_FLAG) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o C62LJ.x
	mv C62LJ.x $(EXE_DIR)/

C62LJ.o : C62LJ.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

C62LJ.f90 : 

clean :
	rm -f *.o
	rm -f *.mod
	rm -f $(EXE_DIR)/C62LJ.x
