#! /bin/bash

source ./config.sh ""
source ./env_vars.sh ""


echo "*** LAPACK ***"
if [ ! -e "$DST_DIR/lib/liblapack.a"  -o "$1" != "" ]; then
	SRC_LINK="https://www.martin-schreiber.info/pub/sweet_local_software/lapack-3.8.0.tar.gz"
	#SRC_LINK="https://www.netlib.org/lapack/lapack-3.6.1.tgz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="lapack-3.8.0"

	cd "$SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1

	tar xzf "$FILENAME"
	cd "$BASENAME"

	ulimit -s 100000 || echo "FAILED TO INCREASE ULIMIT & now ignoring this problem"
	if true; then
		# Use cmake
		mkdir -p build || exit

		cd build || exit
		cmake ../ || exit
		make


		echo "Installing liblapack.a"
		cp ./lib/liblapack.a "$DST_DIR/lib" || exit 1
		#cp ./librefblas.a "$DST_DIR/lib" || exit 1
	else
		# create default configuration
		cp make.inc.example make.inc

		# Replace default Fortran compiler
		if [ ! -z "$FC" ]; then
			# Default C compiler
			sed -i "s/^CC.*=.*/CC = ${CC}/" make.inc

			# Default Fortran compiler
			sed -i "s/^FORTRAN.*=.*/FORTRAN = ${FC}/" make.inc

			# Default linker
			sed -i "s/^LOADER.*=.*/LOADER = ${FC}/" make.inc

			# Get rid of this option
			sed -i "s/-frecursive//" make.inc

			# Use INT_CPU_TIME timer
			sed -i "s/^TIMER.*=.*/TIMER = INT_CPU_TIME/" make.inc

			
		fi

		make blaslib || exit 1
		make || exit 1


		#cd LAPACKE || exit 1
		#make lapacke || exit 1
		#mkdir -p "$DST_DIR/lib"
		echo "Installing liblapack.a"
		cp ./liblapack.a "$DST_DIR/lib" || exit 1
		#cp ./librefblas.a "$DST_DIR/lib" || exit 1
	fi

	echo "DONE"

else
	echo "LAPACK already installed"
fi
