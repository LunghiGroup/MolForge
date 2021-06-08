.SUFFIXES:

SOURCES = cell2abc.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

cell2abc : $(OBJECTS)
	$(FC) $(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_FLAG) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o cell2abc.x
	mv cell2abc.x $(EXE_DIR)/

cell2abc.o : cell2abc.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

cell2abc.f90 : 

clean :
	rm -f *.o
	rm -f *.mod
	rm -f $(EXE_DIR)/cell2abc.x
