#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="SHTNS"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libshtns.a"

PKG_URL_SRC="shtns-3.5.2.tar.gz"


if [ "$CC" == "gcc-6" ]; then
	# Special workaround for gcc-6 compiler
	# This is mainly for the continuous testing environment
	export CFLAGS="$CFLAGS -mno-avx512f"
fi

config_setup

config_package $@


if [ "#$TRAVIS" != "#" ]; then
	echo_info "Detected Travis"
	echo_info "Disabling SIMD because of Travis"
	CONFIGURE_EXTRA_FLAGS="--disable-mkl --disable-knl --disable-cuda --disable-simd"
fi


# Also use special kernel compiler if $CC env variable is set
if [[ ! -z "$CC" ]]; then
	CONFIGURE_EXTRA_FLAGS+=" --enable-kernel-compiler=$CC"
fi


echo_info_hline
echo_info "SHTNS noOpenMP:"
# Python, no OpenMP

config_configure --disable-openmp $CONFIGURE_EXTRA_FLAGS

config_make_clean
config_make_default_install


if [ "`uname`" == "DarwinXXX" ]; then
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


