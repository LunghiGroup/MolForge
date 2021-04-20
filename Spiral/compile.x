### compile all moduli

export EXE_NAME=Spiral.x
export MOLFORGE=/home/Alessandro/ERC/MolForge/src
export MOD_DIR=${MOLFORGE}/Modules
export LIB_DIR=${MOLFORGE}/Modules
export LIBS="-lscalapack -llapack -lmolforge"
export FFLAGS="-fmax-errors=10 -g -O2 -fcheck=all"

mpif90 ${FFLAGS} -c variables.f90 -I ${MOD_DIR}
mpif90 ${FFLAGS} -c parse_input.f90 -I ${MOD_DIR}
mpif90 ${FFLAGS} -c Spiral.f90 -I ${MOD_DIR}
mpif90 ${FFLAGS} *.o -I ${MOD_DIR} -L ${LIB_DIR} ${LIBS} -o ${EXE_NAME}

mv ${EXE_NAME} bin/
