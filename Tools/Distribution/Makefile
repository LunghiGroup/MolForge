.SUFFIXES:

SOURCES = Distribution.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

Distribution : $(OBJECTS)
	$(FC) $(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_FLAG) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o Distribution.x
	mv Distribution.x $(EXE_DIR)/

Distribution.o : Distribution.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

Distribution.f90 : 

clean :
	rm -f *.o
	rm -f *.mod
	rm -f $(EXE_DIR)/Distribution.x
