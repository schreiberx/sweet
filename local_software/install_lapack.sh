#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="lapack"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/liblapack.a"
PKG_URL_SRC="lapack-3.8.0.tar.gz"

config_setup

config_package $@

ulimit -s 100000 || echo_warning "Warning: failed to increase ulimit"

if true; then
	# Use cmake
	mkdir -p build
	cd build

	echo_info "Executing 'cmake'.."
	config_exec cmake ../

	echo_info "Executing 'make'.."
	config_exec make $MAKE_DEFAULT_OPTS


	echo "Installing..."
	config_exec mkdir -p "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib"

	echo " + liblapack.a"
	config_exec cp ./lib/liblapack.a "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib" || echo_error_exit "Failed to install liblapack.a"

	echo " + liblas.a"
	config_exec cp ./lib/libblas.a "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib" || echo_error_exit "Failed to install libblas.a"
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
	
	make $MAKE_DEFAULT_OPTS blaslib
	make $MAKE_DEFAULT_OPTS

	echo "Installing..."
	mkdir -p "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib"

	echo " + liblapack.a"
	cp ./liblapack.a "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib" || echo_error_exit "Failed to install liblapack.a"

	echo " + liblas.a"
	cp ./libblas.a "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib" || echo_error_exit "Failed to install libblas.a"
fi
