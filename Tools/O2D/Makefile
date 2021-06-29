.SUFFIXES:

SOURCES = abc2cell.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

abc2cell : $(OBJECTS)
	$(FC) $(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_FLAG) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o abc2cell.x
	mv abc2cell.x $(EXE_DIR)/

abc2cell.o : abc2cell.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

abc2cell.f90 : 

