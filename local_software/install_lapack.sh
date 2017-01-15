#! /bin/bash

source config.sh


echo "*** LAPACK ***"
if [ ! -e "$DST_DIR/lib/liblapack.a"  -o "$1" != "" ]; then
	SRC_LINK="http://www.martin-schreiber.info/pub/sweet_local_software/lapack-3.6.1.tgz"
	#SRC_LINK="http://www.netlib.org/lapack/lapack-3.6.1.tgz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="lapack-3.6.1"

	cd "$SRC_DIR"

	if [ ! -e "$FILENAME" ]; then
		curl "$SRC_LINK" -o "$FILENAME" || exit 1
	fi

	tar xzf "$FILENAME"
	cd "$BASENAME"

	# create default configuration
	cp make.inc.example make.inc

	#./configure --prefix="$DST_DIR" || exit 1

	ulimit -s 100000 || echo "FAILED TO INCREASE ULIMIT & now ignoring this problem"
	make blaslib || exit 1
	make || exit 1

	#cd LAPACKE || exit 1
	#make lapacke || exit 1
	#mkdir -p "$DST_DIR/lib"
	cp ./liblapack.a "$DST_DIR/lib" || exit 1
	cp ./librefblas.a "$DST_DIR/lib" || exit 1

	echo "DONE"

else
	echo "LAPACK already installed"
fi
