#! /bin/bash

source config.sh


if [ "${HOSTNAME:0:8}" == "cheyenne" ]; then
	echo "Compilation of GCC fails on Cheyenne, please use"
	echo "	module load gnu/7.1.0"
	exit 0
fi

echo "*** GCC7.1 ***"
if [ ! -e "$DST_DIR/bin/gcc-7.1"  -o "$1" != "" ]; then
	SRC_LINK="http://www.martin-schreiber.info/pub/sweet_local_software/gcc-7.1.0.tar.bz2"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="gcc-7.1.0"

	cd "$SRC_DIR"
	download "$SRC_LINK" "$FILENAME" || exit 1

	echo "Uncompressing $FILENAME"
	tar xjf "$FILENAME" || exit 1

	if [ ! -e "$BASENAME" ]; then
		echo "$BASENAME does not exist"
		exit 1
	fi

	##########################
	# GMP
	##########################
	LIB_LINK="https://gmplib.org/download/gmp/gmp-6.1.2.tar.bz2"
	LIB_FILENAME="`basename $LIB_LINK`"
	LIB_BASENAME="gmp-6.1.2"
	LIB_BASENAME_SHORT="gmp"
	download "$LIB_LINK" "$LIB_FILENAME" || exit 1

	if [ ! -e "$BASENAME/$LIB_BASENAME" ]; then
		tar xjf "$LIB_FILENAME" || exit 1
	fi

	mv "$LIB_BASENAME" "$BASENAME/"
	ln -sf "$LIB_BASENAME" "./$BASENAME/$LIB_BASENAME_SHORT"

	##########################
	# MPFR
	##########################
	LIB_LINK="https://ftp.gnu.org/gnu/mpfr/mpfr-3.1.5.tar.bz2"
	LIB_FILENAME="`basename $LIB_LINK`"
	LIB_BASENAME="mpfr-3.1.5"
	LIB_BASENAME_SHORT="mpfr"
	download "$LIB_LINK" "$LIB_FILENAME" || exit 1

	if [ ! -e "$BASENAME/$LIB_BASENAME" ]; then
		tar xjf "$LIB_FILENAME" || exit 1
	fi

	mv "$LIB_BASENAME" "$BASENAME/"
	ln -sf "$LIB_BASENAME" "./$BASENAME/$LIB_BASENAME_SHORT"



	##########################
	# MPC
	##########################
	LIB_LINK="ftp://ftp.gnu.org/gnu/mpc/mpc-1.0.3.tar.gz"
	LIB_FILENAME="`basename $LIB_LINK`"
	LIB_BASENAME="mpc-1.0.3"
	LIB_BASENAME_SHORT="mpc"
	download "$LIB_LINK" "$LIB_FILENAME" || exit 1
	
	if [ ! -e "$BASENAME/$LIB_BASENAME" ]; then
		tar xzf "$LIB_FILENAME" || exit 1
	fi

	mv "$LIB_BASENAME" "$BASENAME/"
	ln -sf "$LIB_BASENAME" "./$BASENAME/$LIB_BASENAME_SHORT"


	##########################
	# ISL
	##########################
	LIB_LINK="http://isl.gforge.inria.fr/isl-0.18.tar.gz"
	LIB_FILENAME="`basename $LIB_LINK`"
	LIB_BASENAME="isl-0.18"
	LIB_BASENAME_SHORT="isl"
	download "$LIB_LINK" "$LIB_FILENAME" || exit 1

	if [ ! -e "$BASENAME/$LIB_BASENAME" ]; then
		tar xzf "$LIB_FILENAME" || exit 1
	fi

	mv "$LIB_BASENAME" "$BASENAME/"
	ln -sf "$LIB_BASENAME" "./$BASENAME/$LIB_BASENAME_SHORT"

	##############################
	##############################

	cd "$BASENAME"
	
	export CC=gcc
	export CXX=g++
	export LINK=ld
	./configure --disable-multilib --enable-languages=c++,fortran  --prefix="$DST_DIR" --program-suffix=-7.1 || exit 1

	make || exit 1
	make install || exit 1

	for i in g++ gcc gcc-ar gcc-nm gcc-ranlib gfortran gcov gcov-tool gfortran; do
		ln -sf "$DST_DIR/bin/$i-7.1" "$DST_DIR/bin/$i"
	done

	echo "DONE"

else
	echo "GCC already installed"
fi
