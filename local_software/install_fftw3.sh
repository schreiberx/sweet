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

config_package $@ || exit 1

CONF_FLAGS=""

if [ "`uname -s`" == "Linux" -o "`uname -s`" == "Darwin" ]; then
	CONF_FLAGS+=" --enable-openmp "
fi

# Activate vectorization code

#	sse only works with single precision
#	CONF_FLAGS+=" --enable-sse"

#
# The runtime autodetect feature in FFTW seems to be buggy
#
# Therefore we use our own autodetection
#
echo_info_hline
echo_info "Autodetect CPU features"
echo_info_hline
CPUFLAGS=$(cat /proc/cpuinfo  | grep ^flags | head -n 1 | sed "s/.*: //")
FEATURES="sse2 avx avx2 avx512"
for FEATURE in $FEATURES; do
	if [[ $CPUFLAGS =~ .*$FEATURE.* ]]; then
		echo_info "Detected '${FEATURE}' in CPU flags"
		CONF_FLAGS+=" --enable-$FEATURE"
	else
		echo_warning "Feature '${FEATURE}' not found in CPU flags"
	fi
done
echo_info_hline

# Never used directly in Fortran
CONF_FLAGS+=" --disable-fortran"

# Enable generation of shared library
CONF_FLAGS+=" --enable-shared"

#	neon only works with single precision
#	CONF_FLAGS+=" --enable-neon"

echo "Configuration flags: $CONF_FLAGS"

config_configure $CONF_FLAGS || exit 1

config_make_default_install || exit 1

config_success || exit 1
