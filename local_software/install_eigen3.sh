#! /bin/bash

source config.sh


echo "*** EIGEN3 ***"
if [ ! -e "$DST_DIR/include/librdic.so.0"  -o "$1" != "" ]; then
	SRC_LINK="http://www.martin-schreiber.info/pub/sweet_local_software/eigen-3.3.0.tar.bz2"
	#SRC_LINK="https://bitbucket.org/eigen/eigen/get/3.3.0.tar.bz2"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="eigen-3.3.0"

	cd "$SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1
	tar xjf "$FILENAME"
	mv "eigen-eigen-26667be4f70b" "$BASENAME"
	cd "$BASENAME"

	mkdir -p build
	cd build
	cmake ../  -DCMAKE_INSTALL_PREFIX="$DST_DIR" || exit 1
	make install || exit 1

	echo "DONE"

else
	echo "RDIC is already installed"
fi
