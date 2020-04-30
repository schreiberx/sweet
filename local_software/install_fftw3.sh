#! /bin/bash

source ./install_helpers.sh ""
source ./env_vars.sh ""

PKG_NAME="fftw"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libfftw3.a"
PKG_URL_SRC="https://www.martin-schreiber.info/pub/sweet_local_software/fftw-3.3.8.tar.gz"

config_setup

config_package $@


CONF_FLAGS=""

if [ "`uname -s`" == "Linux" -o "`uname -s`" == "Darwin" ]; then
	CONF_FLAGS+=" --enable-openmp "
fi

# Activate vectorization code

#	sse only works with single precision
#	CONF_FLAGS+=" --enable-sse"

if [ "`uname -s`" == "x86_64" ]; then
	CONF_FLAGS+=" --enable-sse2"
	CONF_FLAGS+=" --enable-avx"
	CONF_FLAGS+=" --enable-avx2"
	CONF_FLAGS+=" --enable-avx512"
	#CONF_FLAGS+=" --enable-avx-128-fma"
fi

# Never used directly in Fortran
#CONF_FLAGS+=" --disable-fortran"
CONF_FLAGS+=" --enable-fortran"

# Enable generation of shared library
CONF_FLAGS+=" --enable-shared"

#	neon only works with single precision
#	CONF_FLAGS+=" --enable-neon"

echo "Configuration flags: $CONF_FLAGS"

config_configure $CONF_FLAGS

config_make_default_install

config_success
