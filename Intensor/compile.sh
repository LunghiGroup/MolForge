### compile all moduli

export EXE_NAME=Intensor.x
export MOLFORGE=/home/Alessandro/ERC/MolForge
export MOD_DIR=${MOLFORGE}/Modules
export LIB_DIR=${MOLFORGE}/libs
export LIBS="-lMolForge -lscalapack -llapack"
export FFLAGS="-fmax-errors=10 -g -O2 -fcheck=all"

mpif90 ${FFLAGS} -c Intensor.f90 -I ${MOD_DIR}
mpif90 ${FFLAGS} *.o -I ${MOD_DIR} -L ${LIB_DIR} ${LIBS} -o ${EXE_NAME}

mv ${EXE_NAME} bin/
