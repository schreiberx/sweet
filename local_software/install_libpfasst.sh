#! /bin/bash

source config.sh


echo "*** libPFASST ***"
if [ ! -e "$DST_DIR/lib/libpfasst.a"  -o "$1" != "" ]; then
	#BASENAME="libpfasst"
	BASENAME="libpfasst_martin"

	cd "$SRC_DIR"

	tar xzf libpfasst_2018_04_23_francois.tar.gz

	#if [ ! -e "$BASENAME" ]; then
	#	git clone https://bitbucket.org/memmett/libpfasst.git || exit 1
	#fi

	cd "$BASENAME"

	sed -i "s/ftn/mpif90/" Makefile.defaults || exit 1

	make clean || exit 1
	make || exit 1

	mkdir -p "$DST_DIR/lib/"

	echo cp -v -f "./lib/libpfasst.a" "$DST_DIR/lib/"
	cp -v -f "./lib/libpfasst.a" "$DST_DIR/lib/"

	echo cp -v -f ./include/* "$DST_DIR/include/"
	cp -v -f ./include/* "$DST_DIR/include/"

	echo "DONE"

else
	echo "libpfasst (Fortran PFASST) is already installed"
fi
