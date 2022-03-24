.SUFFIXES:

SOURCES = O2D.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

O2D : $(OBJECTS)
	$(FC) $(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_FLAG) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o O2D.x
	mv O2D.x $(EXE_DIR)/

O2D.o : O2D.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

O2D.f90 : 

clean :
	rm -f *.o
	rm -f *.mod
	rm -f $(EXE_DIR)/O2D.x
