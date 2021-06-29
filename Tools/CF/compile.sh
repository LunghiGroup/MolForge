	export MOLFORGE=/home/Alessandro/ERC/MolForge/src
	export MOD_DIR=${MOLFORGE}/Modules
	export LIB_DIR=${MOLFORGE}/Modules
	mpif90 -g -O2 *f90 -o CF.x -I ${MOD_DIR} -L ${LIB_DIR} -lscalapack -llapack -lmolforge -fmax-errors=10
	
