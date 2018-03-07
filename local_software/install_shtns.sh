#! /bin/bash

source config.sh


echo "*** SHTNS ***"
if [ ! -e "$DST_DIR/lib/libshtns_omp.a"  -o "$1" != "" ]; then
	SRC_LINK="https://www.martin-schreiber.info/pub/sweet_local_software/shtns-3.0-r618.tar.gz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="shtns-3.0-r618"

	cd "$SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1

	tar xzf "$FILENAME"
	cd "$BASENAME"

	# no OpenMP
	./configure --prefix="$DST_DIR" --enable-python --disable-openmp || exit 1
	make || exit 1
	#make install || exit 1
	python3 setup.py install --prefix="$DST_DIR"

	# with OpenMP
	#make clean
	./configure --prefix="$DST_DIR" --enable-python --enable-openmp || exit 1
	make || exit 1
	python3 setup.py install --prefix="$DST_DIR"
	#make install || exit 1

	echo "DONE"

else
	echo "SHTNS already installed"
fi
