	export MOLFORGE=$( pwd )
	export SRC_DIR=${MOLFORGE}/src
	export MOD_DIR=${SRC_DIR}/Modules
	export LIB_DIR=${SRC_DIR}/Modules
	cd ${MOD_DIR}
###        ./compile.x
	cd ${SRC_DIR}
	mpif90 -g -O3 *f90 -o MolForge.x -I ${MOD_DIR} -L ${LIB_DIR} -lscalapack -llapack -lmolforge -fmax-errors=10
	cd ${MOLFORGE}


	
