#! /bin/bash

source ./install_helpers.sh "" || exit 1
source ./env_vars.sh "" || exit 1

# Name of package
PKG_NAME="fftw"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libfftw3.a"

# URL to source code to fetch it
PKG_URL_SRC="fftw-3.3.8.tar.gz"

# subdirectory of source in extracted package
# (autodetected with basename of url without file extension if not set)
#PKG_SRC_SUBDIR=""

config_package $@

CONF_FLAGS=""

if [ "`uname -s`" == "Linux" -o "`uname -s`" == "Darwin" ]; then
	CONF_FLAGS+=" --enable-openmp "
fi

# Activate vectorization code

#	sse only works with single precision
#	CONF_FLAGS+=" --enable-sse"

CONF_FLAGS+=" --enable-sse2"
CONF_FLAGS+=" --enable-avx"
CONF_FLAGS+=" --enable-avx2"
CONF_FLAGS+=" --enable-avx512"
#CONF_FLAGS+=" --enable-avx-128-fma"

# Never used directly in Fortran
CONF_FLAGS+=" --disable-fortran"

# Enable generation of shared library
CONF_FLAGS+=" --enable-shared"

#	neon only works with single precision
#	CONF_FLAGS+=" --enable-neon"

echo "Configuration flags: $CONF_FLAGS"

config_configure $CONF_FLAGS

config_make_default_install

config_success
