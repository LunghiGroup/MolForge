 	export LIB_DIR=/home/Alessandro/ERC/MolForge/src/Modules/ 
 	export INCLUDE_DIR=/home/Alessandro/ERC/MolForge/src/Modules/ 

	gfortran rotmat.f90 -I ${INCLUDE_DIR} -L ${LIB_DIR} -llapack -lmolforge -o rotmat.x
