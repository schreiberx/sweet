#! /bin/bash

source config.sh ""
source env_vars.sh ""


echo "*** SHTNS ***"
SRC_LINK="https://www.martin-schreiber.info/pub/sweet_local_software/nschaeff-shtns-2018_10_01.tar.bz2"
FILENAME="`basename $SRC_LINK`"
BASENAME="nschaeff-shtns-2018_10_01"

if [ ! -e "$DST_DIR/lib/libshtns_omp.a"  -o "$1" != "" ]; then

	cd "$SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1

	tar xjf "$FILENAME"
	cd "$BASENAME"

	
	# library, OpenMP
	make clean
	./configure --prefix="$DST_DIR" --enable-openmp || exit 1
	make install || exit 1

	# library, no OpenMP
	make clean
	./configure --prefix="$DST_DIR" --disable-openmp || exit 1
	make install || exit 1


	echo "DONE"

else
	echo "SHTNS already installed"
fi
