#! /bin/bash

source ./config.sh ""
source ./env_vars.sh ""


echo "*** libPFASST ***"

if [ ! -e "$DST_DIR/lib/libpfasst.a"  -o "$1" != "" ]; then
	SRC_LINK="https://www.martin-schreiber.info/pub/sweet_local_software/libpfasst_sweet_2018_09_27.tar.bz2"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="libpfasst_sweet_2018_09_27"

	cd "$SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1
	tar xjf "$FILENAME"
	cd "$BASENAME"
	sed -i "s/ftn/mpif90/" Makefile.defaults || exit 1

	make clean || exit 1
	make FC=$SWEET_MPIF90 || exit 1

	mkdir -p "$DST_DIR/lib/"

	# Copy static library
	echo cp -v -f "./lib/libpfasst.a" "$DST_DIR/lib/"
	cp -v -f "./lib/libpfasst.a" "$DST_DIR/lib/"

	# Copy modules
	echo cp -v -f ./include/*mod "$DST_DIR/include/"
	cp -v -f ./include/*mod "$DST_DIR/include/"

	echo "DONE"

else
	echo "libpfasst (Fortran PFASST) is already installed"
fi
