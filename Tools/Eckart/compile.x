	export MOLFORGE=/home/Alessandro/ERC/MolForge/src
	export MOD_DIR=${MOLFORGE}/Modules
	export LIB_DIR=${MOLFORGE}/Modules
	mpif90 -g -o Eckart.x -I ${MOD_DIR} -L ${LIB_DIR} Eckart.f90 -llapack -lmolforge 
