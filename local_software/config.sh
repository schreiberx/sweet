
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

function download {
	echo "Downloading from $1"
	#curl -C - "$1" -o "$2" || exit 1
	wget --continue --progress=bar "$1" -O "$2" || exit 1
}
