#! /bin/bash

source ./config_install.sh ""
source ./env_vars.sh ""


echo "*** TBB2 ***"
if [ -e "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libtbb.so"  -o "$1" != "" ]; then
	echo "TBB2 already installed"
	exit 1
fi

SRC_LINK="https://www.martin-schreiber.info/pub/sweet_local_software/tbb2018_20170919oss_lin.tgz"
#SRC_LINK="https://www.fftw.org/fftw-3.3.4.tar.gz"
FILENAME="`basename $SRC_LINK`"
#	BASENAME="fftw-3.3.4"
BASENAME="tbb2018_20170919oss"

cd "$SWEET_LOCAL_SOFTWARE_SRC_DIR"

download "$SRC_LINK" "$FILENAME" || exit 1
tar xzf "$FILENAME"
cd "$BASENAME"


cp -r ./bin ./include ./python ./lib "$SWEET_LOCAL_SOFTWARE_DST_DIR"
#sed -i "s/SUBSTITUTE_INSTALL_DIR_HERE/auto_tbbroot/" "$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/tbbvars.sh"

echo "DONE"
