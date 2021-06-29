.SUFFIXES:

SOURCES = Mat2Euler.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

Mat2Euler : $(OBJECTS)
	$(FC) $(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_FLAG) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o Mat2Euler.x
	mv Mat2Euler.x $(EXE_DIR)/

Mat2Euler.o : Mat2Euler.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

Mat2Euler.f90 : 

clean :
	rm -f *.o
	rm -f *.mod
	rm -f $(EXE_DIR)/Mat2Euler.x
