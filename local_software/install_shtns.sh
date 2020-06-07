#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="SHTNS"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libshtns.a"

# This version doesn't work on Travis CI
PKG_URL_SRC="shtns-3.3.1-r694.tar.gz"
#PKG_URL_SRC="shtns-3.1-r633.tar.gz"

config_setup

config_package $@

echo_info_hline
echo_info "SHTNS noOpenMP:"
# Python, no OpenMP

if [ "#$TRAVIS" != "#" ]; then
	echo_info "Detected Travis"
	echo_info "Disabling SIMD because of Travis"
	EXTRA_FLAGS="--disable-mkl --disable-knl --disable-cuda --disable-simd"
fi

# Special flag for sk2 (@ CAPS hardware)
if [ "#$(hostname)" = "#sk1" -o "#$(hostname)" = "#sk2" ]; then
       export CFLAGS="$CFLAGS -march=skylake"
fi

config_configure --disable-openmp $EXTRA_FLAGS

# Special flag for sk2 (@ CAPS hardware)
pwd
if [ "#$(hostname)" = "#sk1" -o "#$(hostname)" = "#sk2" ]; then
	sed -i "s/-march=native/-march=skylake/" "Makefile"
fi

config_make_clean
config_make_default_install

echo_info_hline
echo_info "SHTNS OpenMP:"
# Python, OpenMP
config_configure --enable-openmp $EXTRA_FLAGS

# Special flag for sk2 (@ CAPS hardware)
if [ "#$(hostname)" = "#sk1" -o "#$(hostname)" = "#sk2" ]; then
	sed -i "s/-march=native/-march=skylake/" "Makefile"
fi

config_make_clean
config_make_default_install

config_success
