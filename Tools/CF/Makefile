.SUFFIXES:

SOURCES = CF.f90

OBJECTS = $(subst .f90,.o,$(SOURCES))

CF : $(OBJECTS) get_SH_Orca get_CF_Orca get_DCF_Orca get_DSH_Orca
	$(FC) $(OBJECTS) $(MOLFORGE_INC_FLAG) $(INC_FLAG) $(MOLFORGE_LINK_FLAG) $(LINK_FLAG) -o CF.x
	mv CF.x $(EXE_DIR)/

CF.o : CF.f90
	$(FC) $(FFLAGS) -c $(subst .o,.f90,$@) $(MOLFORGE_INC_FLAG) $(INC_FLAG)

CF.f90 : 

get_CF_Orca :
	sed "s/XXXSETDIRXXX/$(subst /,\/,$(MOLFORGE_DIR))/g" get_CF_Orca.x > $(EXE_DIR)/get_CF_Orca.x
	chmod 777 $(EXE_DIR)/get_CF_Orca.x

get_DCF_Orca :
	sed "s/XXXSETDIRXXX/$(subst /,\/,$(MOLFORGE_DIR))/g" get_DCF_Orca.x > $(EXE_DIR)/get_DCF_Orca.x
	chmod 777 $(EXE_DIR)/get_DCF_Orca.x

get_SH_Orca :
	sed "s/XXXSETDIRXXX/$(subst /,\/,$(MOLFORGE_DIR))/g" get_SH_Orca.x > $(EXE_DIR)/get_SH_Orca.x
	chmod 777 $(EXE_DIR)/get_SH_Orca.x

get_DSH_Orca :
	sed "s/XXXSETDIRXXX/$(subst /,\/,$(MOLFORGE_DIR))/g" get_DSH_Orca.x > $(EXE_DIR)/get_DSH_Orca.x
	chmod 777 $(EXE_DIR)/get_DSH_Orca.x

clean :
	rm -f *.o
	rm -f *.mod
	rm -f $(EXE_DIR)/CF.x
	rm -f $(EXE_DIR)/get_SH_Orca.x
	rm -f $(EXE_DIR)/get_CF_Orca.x
	rm -f $(EXE_DIR)/get_DSH_Orca.x
	rm -f $(EXE_DIR)/get_DCF_Orca.x
