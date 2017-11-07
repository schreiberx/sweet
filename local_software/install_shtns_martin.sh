#! /bin/bash

source config.sh


echo "*** SHTNS ***"
if [ ! -e "$DST_DIR/lib/libshtns_omp.a"  -o "$1" != "" ]; then
	#SRC_LINK="http://www.martin-schreiber.info/pub/sweet_local_software/shtns-2.6.6-r535.tar.gz"
	#SRC_LINK="https://bitbucket.org/nschaeff/shtns/downloads/shtns-2.6.6-r535.tar.gz"
	#FILENAME="`basename $SRC_LINK`"
	FILENAME="shtns_martin_2017_11_06.tar.gz"
	BASENAME="shtns_martin"

	cd "$SRC_DIR"

	#download "$SRC_LINK" "$FILENAME" || exit 1

	tar xzf "$FILENAME"
	cd "$BASENAME"

	# no OpenMP
	make clean
	./configure --prefix="$DST_DIR" --disable-openmp || exit 1
	make install || exit 1

	# with OpenMP
	make clean
	./configure --prefix="$DST_DIR" --enable-openmp || exit 1
	make install || exit 1

	echo "DONE"

else
	echo "SHTNS already installed"
fi
