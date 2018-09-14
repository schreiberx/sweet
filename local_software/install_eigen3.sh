#! /bin/bash

source config.sh
source env_vars.sh


echo "*** EIGEN3 ***"
if [ ! -e "$DST_DIR/include/eigen3/Eigen/Eigenvalues"  -o "$1" != "" ]; then
	SRC_LINK="https://www.martin-schreiber.info/pub/sweet_local_software/eigen-3.3.3.tar.bz2"
	#SRC_LINK="https://bitbucket.org/eigen/eigen/get/3.3.3.tar.bz2"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="eigen-3.3.3"

	cd "$SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1
	tar xjf "$FILENAME"
	mv "eigen-eigen-67e894c6cd8f" "$BASENAME"
	cd "$BASENAME"

	mkdir -p build
	cd build
	cmake ../  -DCMAKE_INSTALL_PREFIX="$DST_DIR" || exit 1
	make install || exit 1

	echo "DONE"

else
	echo "Eigen is already installed"
fi
