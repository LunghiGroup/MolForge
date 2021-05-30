### compile all moduli

export EXE_NAME=Spiral.x
export SPIRAL_DIR=$(pwd) 
export MOLFORGE=$( echo ${SPIRAL_DIR} | sed "s/Spiral/src/g" )
export MOD_DIR=${MOLFORGE}/Modules
export LIB_DIR=${MOLFORGE}/Modules

cd ${MOD_DIR}
./compile.x
cd ${SPIRAL_DIR}

export LIBS="-lscalapack -llapack -lmolforge"
export FFLAGS="-fmax-errors=10 -g -O2 -fcheck=all"

rm *.o *.mod

mpif90 ${FFLAGS} -c variables.f90 -I ${MOD_DIR}
mpif90 ${FFLAGS} -c parse_input.f90 -I ${MOD_DIR}
mpif90 ${FFLAGS} -c Spiral.f90 -I ${MOD_DIR}
mpif90 ${FFLAGS} *.o -I ${MOD_DIR} -L ${LIB_DIR} ${LIBS} -o ${EXE_NAME}

mv ${EXE_NAME} bin/
