export MODULEFILE=$PWD/../ModelBuild/ogstm/compilers/machine_modules/g100.intel
source $MODULEFILE



OGSTM_ARCH=x86_64
OGSTM_OS=LINUX
OGSTM_COMPILER=intel


DEBUG_OCEANVAR=.dbg


INC_FILE=${OGSTM_ARCH}.${OGSTM_OS}.${OGSTM_COMPILER}${DEBUG_OCEANVAR}.inc
cp $INC_FILE compiler.inc
make clean
gmake
