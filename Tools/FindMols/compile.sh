	export MOLFORGE=/home/Alessandro/ERC/MolForge/src
	export MOD_DIR=${MOLFORGE}/Modules
	export LIB_DIR=${MOLFORGE}/Modules
	mpif90 -g -O2 *f90 -o find_mols.x -I ${MOD_DIR} -L ${LIB_DIR} -lmolforge -llapack -lscalapack -fmax-errors=10
	
