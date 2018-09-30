#! /bin/bash

source ./config.sh ""
source ./env_vars.sh ""


echo "*** OpenMPI ***"
if [ ! -e "$DST_DIR/bin/mpicc"  -o "$1" != "" ]; then
	SRC_LINK="https://www.martin-schreiber.info/pub/sweet_local_software/openmpi-1.10.2.tar.bz2"
	#SRC_LINK="https://www.open-mpi.org/software/ompi/v1.10/downloads/openmpi-1.10.2.tar.bz2"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="openmpi-1.10.2"

	cd "$SRC_DIR"
	download "$SRC_LINK" "$FILENAME" || exit 1

	echo "Uncompressing $FILENAME"
	tar xjf "$FILENAME"

	if [ ! -e "$BASENAME" ]; then
		echo "$BASENAME does not exist"
		exit 1
	fi


	##############################
	##############################

	cd "$BASENAME"
	
	export CC=gcc
	export CXX=g++
	export LINK=ld
	./configure --enable-mpi-fortran --prefix="$DST_DIR" || exit 1

	make || exit 1
	make install || exit 1

	echo "DONE"

else
	echo "OpenMPI already installed"
fi
