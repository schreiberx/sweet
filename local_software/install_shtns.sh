#! /bin/bash

source config.sh


echo "*** SHTNS ***"
if [ ! -e "$DST_DIR/lib/libshtns_omp.a"  -o "$1" != "" ]; then
	SRC_LINK="https://www.martin-schreiber.info/pub/sweet_local_software/shtns-2.9_alpha.tar.gz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="shtns-2.9_alpha"

	cd "$SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1

	tar xzf "$FILENAME"
	cd "$BASENAME"

	# no OpenMP
	make clean
	./configure --prefix="$DST_DIR" --enable-python --disable-openmp || exit 1
	python setup.py install --prefix="$DST_DIR"
	make install || exit 1

	# with OpenMP
	make clean
	./configure --prefix="$DST_DIR" --enable-python --enable-openmp || exit 1
	python setup.py install --prefix="$DST_DIR"
	make install || exit 1

	echo "DONE"

else
	echo "SHTNS already installed"
fi
