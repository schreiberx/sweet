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
	EXTRA_FLAGS="--disable-mkl --disable-knl --disable-cuda --disable-avx512"
fi

config_configure --disable-mem --disable-openmp $EXTRA_FLAGS
config_make_clean
config_make_default_install

echo_info_hline
echo_info "SHTNS OpenMP:"
# Python, OpenMP
config_configure --disable-mem --enable-openmp $EXTRA_FLAGS
#config_configure --enable-openmp
config_make_clean
config_make_default_install

config_success
