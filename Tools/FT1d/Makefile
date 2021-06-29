.SUFFIXES:

SOURCES = FT.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

FT1d : $(OBJECTS)
	$(FC) $(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_FLAG) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o FT1d.x
	mv FT1d.x $(EXE_DIR)/

FT.o : FT.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

FT.f90 : 

clean :
	rm -f *.o
	rm -f *.mod
	rm -f $(EXE_DIR)/FT1d.x
