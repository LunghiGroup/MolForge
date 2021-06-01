### compile all moduli

export EXE_NAME=PhonDy.x
export PHONDY_DIR=$(pwd) 
export MOLFORGE=$( echo ${PHONDY_DIR} | sed "s/PhonDy/src/g" )
export MOD_DIR=${MOLFORGE}/Modules
export LIB_DIR=${MOLFORGE}/Modules

cd ${MOD_DIR}
./compile.x
cd ${PHONDY_DIR}

export LIBS="-lscalapack -llapack -lmolforge"
export FFLAGS="-fmax-errors=10 -g -O2 -fcheck=all"

mpif90 ${FFLAGS} -c PhonDy.f90 -I ${MOD_DIR}
mpif90 ${FFLAGS} *.o -I ${MOD_DIR} -L ${LIB_DIR} ${LIBS} -o ${EXE_NAME}

mkdir -p bin/
mv ${EXE_NAME} bin/
