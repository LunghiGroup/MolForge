.SUFFIXES:

SOURCES = get_bispectrum.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

Bispectrums : $(OBJECTS)
	$(FC) $(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_FLAG) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o get_bis.x
	mv get_bis.x $(EXE_DIR)/

get_bispectrum.o : get_bispectrum.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

get_bispectrum.f90 : 

clean :
	rm -f *.o
	rm -f *.mod
	rm -f $(EXE_DIR)/get_bis.x
