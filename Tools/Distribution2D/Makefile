.SUFFIXES:

SOURCES = Distribution2D.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

Distribution2D : $(OBJECTS)
	$(FC) $(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_FLAG) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o Distribution2D.x
	mv Distribution2D.x $(EXE_DIR)/

Distribution2D.o : Distribution2D.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

Distribution2D.f90 : 

clean :
	rm -f *.o
	rm -f *.mod
	rm -f $(EXE_DIR)/Distribution2D.x
