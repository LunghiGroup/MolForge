.SUFFIXES:

export INC_DIR = 
export LIB_DIR = 
export LIBS = -lscalapack -llapack 

export FC = mpif90
export FFLAGS = -g -O2
export INC_FLAG =  ## -I$(INC_DIR) 
export LINK_FLAG = $(LIBS) ## -L$(LIB_DIR) $(LIBS)

export SCALAPACK_FLAG = GNU ## Has to be "MKL" or "GNU", any other value will default to GNU


######## The rest of the flags should not be touched

export MOLFORGE_DIR = $(shell pwd )

export MODULES_DIR = $(MOLFORGE_DIR)/Modules

export SPIRAL_DIR = $(MOLFORGE_DIR)/Spiral
export PHONDY_DIR = $(MOLFORGE_DIR)/PhonDy
export TOOLS_DIR = $(MOLFORGE_DIR)/Tools

export EXE_DIR = $(MOLFORGE_DIR)/bin

export MOLFORGE_INC_DIR = $(MODULES_DIR)
export MOLFORGE_LIB_DIR = $(MOLFORGE_DIR)/libs
export MOLFORGE_LIBS = -lMolForge 

export MOLFORGE_INC_FLAG = -I$(MOLFORGE_INC_DIR) 
export MOLFORGE_LINK_FLAG = -L$(MOLFORGE_LIB_DIR) $(MOLFORGE_LIBS)


all : libMolForge Spiral PhonDy Tools


libMolForge :
	mkdir -p $(MOLFORGE_DIR)/libs
	cd $(MODULES_DIR) && $(MAKE) $@

Spiral : libMolForge
	mkdir -p $(MOLFORGE_DIR)/bin
	cd $(SPIRAL_DIR) && $(MAKE) $@

PhonDy : libMolForge
	mkdir -p $(MOLFORGE_DIR)/bin
	cd $(PHONDY_DIR) && $(MAKE) $@

Tools : libMolForge
	mkdir -p $(MOLFORGE_DIR)/bin
	cd $(TOOLS_DIR) && $(MAKE) all


clean : clean_MolForge clean_Spiral clean_PhonDy clean_Tools

clean_MolForge :
	rm -f $(MODULES_DIR)/*.o
	rm -f $(MODULES_DIR)/*.mod
	rm -f $(MOLFORGE_LIB_DIR)/libMolForge.a

clean_Spiral : 
	rm -f $(SPIRAL_DIR)/*.o
	rm -f $(SPIRAL_DIR)/*.mod
	rm -f $(EXE_DIR)/Spiral.x

clean_PhonDy : 
	rm -f $(PHONDY_DIR)/*.o
	rm -f $(PHONDY_DIR)/*.mod
	rm -f $(EXE_DIR)/PhonDy.x

clean_Tools :
	cd $(TOOLS_DIR) && $(MAKE) clean

