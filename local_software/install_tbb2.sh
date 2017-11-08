#! /usr/bin/env bash

source config.sh


echo "*** TBB2 ***"
if [ -e "$DST_DIR/lib/libtbb.so"  -o "$1" != "" ]; then
	echo "TBB2 already installed"
	exit 1
fi

SRC_LINK="http://www.martin-schreiber.info/pub/sweet_local_software/tbb2018_20170919oss_lin.tgz"
#SRC_LINK="http://www.fftw.org/fftw-3.3.4.tar.gz"
FILENAME="`basename $SRC_LINK`"
#	BASENAME="fftw-3.3.4"
BASENAME="tbb2018_20170919oss"

cd "$SRC_DIR"

download "$SRC_LINK" "$FILENAME" || exit 1
tar xzf "$FILENAME"
cd "$BASENAME"


cp -r ./bin ./include ./python ./lib "$DST_DIR"
#sed -i "s/SUBSTITUTE_INSTALL_DIR_HERE/auto_tbbroot/" "$DST_DIR/bin/tbbvars.sh"

echo "DONE"
