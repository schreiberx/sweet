
# dont touch this!

PWD=`pwd`
ROOT_DIR="$PWD/local"
SRC_DIR="$PWD/local_src"
DST_DIR="$ROOT_DIR"

mkdir -p "$ROOT_DIR"
mkdir -p "$SRC_DIR"

# Intel compiler flags
#export CC=icc
#export CXX=icpc
#export LINK=icpc


# E.g. for PowerPC
#export CC=xlc
#export CXX=xlc++
#export LINK=xlc++

# on BGQ to compile FFT
#export HOST="powerpc64-bgq-linux"

