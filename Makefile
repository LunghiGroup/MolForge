.SUFFIXES:

export MOLFORGE_DIR = $(shell pwd )
export MODULES_DIR = $(MOLFORGE_DIR)/Modules
export SPIRAL_DIR = $(MOLFORGE_DIR)/Spiral
export PHONDY_DIR = $(MOLFORGE_DIR)/PhonDy
export TOOLS_DIR = $(MOLFORGE_DIR)/Tools

export INC_DIR = $(MODULES_DIR)
export LIB_DIR = $(MODULES_DIR)
export LIBS = -lmolforge -lscalapack -llapack -lblas

export FC = mpif90
export FFLAGS = -g -fmax-errors=10 -O2
export FFLAGS += -I$(INC_DIR) -L$(LIB_DIR) $(LIBS)

export SCALAPACK_FLAG = GNU ## Has to be "MKL" or "GNU", any other value will default to GNU

libMolForge :
	cd $(MODULES_DIR) && $(MAKE) libMolForge

##Spiral : libMolForge
##	cd $(SPIRAL_DIR) && $(MAKE) 

##PhonDy : libMolForge
##	cd $(PHONDY_DIR) && $(MAKE) 

