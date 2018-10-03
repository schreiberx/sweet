#! /bin/bash

source ./config_install.sh ""
source ./env_vars.sh ""


echo "*** libPFASST ***"

if [ ! -e "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libpfasst.a"  -o "$1" != "" ]; then
	SRC_LINK="https://www.martin-schreiber.info/pub/sweet_local_software/libpfasst_sweet_2018_09_27.tar.bz2"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="libpfasst_sweet_2018_09_27"

	cd "$SWEET_LOCAL_SOFTWARE_SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1
	tar xjf "$FILENAME"
	cd "$BASENAME"
	sed -i "s/ftn/mpif90/" Makefile.defaults || exit 1

	make clean || exit 1
	make FC=$SWEET_MPIF90 || exit 1

	mkdir -p "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/"

	# Copy static library
	echo cp -v -f "./lib/libpfasst.a" "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/"
	cp -v -f "./lib/libpfasst.a" "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/"

	# Copy modules
	echo cp -v -f ./include/*mod "$SWEET_LOCAL_SOFTWARE_DST_DIR/include/"
	cp -v -f ./include/*mod "$SWEET_LOCAL_SOFTWARE_DST_DIR/include/"

	echo "DONE"

else
	echo "libpfasst (Fortran PFASST) is already installed"
fi
