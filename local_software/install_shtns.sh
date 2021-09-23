#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="SHTNS"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libshtns.a"

PKG_URL_SRC="shtns-3.4.6.tar.gz"

config_setup

config_package $@

if [ "#$TRAVIS" != "#" ]; then
	echo_info "Detected Travis"
	echo_info "Disabling SIMD because of Travis"
	CONFIGURE_EXTRA_FLAGS="--disable-mkl --disable-knl --disable-cuda --disable-simd"
fi

#CONFIGURE_EXTRA_FLAGS+=" --enable-ishioka"
#CONFIGURE_EXTRA_FLAGS+=" --disable-ishioka"


echo_info_hline
echo_info "SHTNS noOpenMP:"
# Python, no OpenMP

config_configure --disable-openmp $CONFIGURE_EXTRA_FLAGS

config_make_clean
config_make_default_install


if [ "`uname`" == "Darwin" ]; then
	echo_info_hline
	echo_info "SHTNS Skipping OpenMP due to lack of OpenMP compiler with clang"
	echo_info_hline

	config_success
	return
fi


echo_info_hline
echo_info "SHTNS OpenMP:"
# Python, OpenMP
config_configure --enable-openmp $CONFIGURE_EXTRA_FLAGS


config_make_clean
config_make_default_install

config_success


