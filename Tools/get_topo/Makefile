.SUFFIXES:

SOURCES = pbc.f90 get_topo.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

get_topo : $(OBJECTS)
	$(FC) $(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_FLAG) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o get_topo.x
	mv get_topo.x $(EXE_DIR)/

get_topo.o : pbc.o get_topo.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

get_topo.f90 : 

pbc.o : pbc.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

pbc.f90 : 

clean :
	rm -f *.o
	rm -f *.mod
	rm -f $(EXE_DIR)/get_topo.x
