.SUFFIXES:

export MOLFORGE_DIR = $(shell pwd )
export MODULES_DIR = $(MOLFORGE_DIR)/Modules
export SPIRAL_DIR = $(MOLFORGE_DIR)/Spiral
export PHONDY_DIR = $(MOLFORGE_DIR)/PhonDy
export TOOLS_DIR = $(MOLFORGE_DIR)/Tools
export EXE_DIR = $(MOLFORGE_DIR)/bin

export INC_DIR =
export LIB_DIR =
export LIBS = -lscalapack -llapack 

export MOLFORGE_INC_DIR = $(MOLFORGE_DIR)/libs
export MOLFORGE_LIB_DIR = $(MOLFORGE_DIR)/libs
export MOLFORGE_LIBS = -lMolForge 

export FC = mpif90
export FFLAGS = -g -fmax-errors=10 -O2
export FFLAGS += -I$(INC_DIR) -L$(LIB_DIR) $(LIBS)

export SCALAPACK_FLAG = GNU ## Has to be "MKL" or "GNU", any other value will default to GNU

libMolForge :
	cd $(MODULES_DIR) && $(MAKE) $@

Spiral : libMolForge
	cd $(SPIRAL_DIR) && $(MAKE) $@

##PhonDy : libMolForge
##	cd $(PHONDY_DIR) && $(MAKE) 

