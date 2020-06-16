#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="SHTNS"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libshtns.a"

PKG_URL_SRC="shtns-3.4.tar.gz"

config_setup

config_package $@

if [ "#$TRAVIS" != "#" ]; then
	echo_info "Detected Travis"
	echo_info "Disabling SIMD because of Travis"
	CONFIGURE_EXTRA_FLAGS="--disable-mkl --disable-knl --disable-cuda --disable-simd"
fi

CONFIGURE_EXTRA_FLAGS+=" --enable-ishioka"


echo_info_hline
echo_info "SHTNS noOpenMP:"
# Python, no OpenMP

# Special flag for sk2 (@ CAPS hardware)
if [ "#$(hostname)" = "#sk1" -o "#$(hostname)" = "#sk2" ]; then
       export CFLAGS="$CFLAGS -march=skylake"
fi

config_configure --disable-openmp $CONFIGURE_EXTRA_FLAGS

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
config_configure --enable-openmp $CONFIGURE_EXTRA_FLAGS

# Special flag for sk2 (@ CAPS hardware)
if [ "#$(hostname)" = "#sk1" -o "#$(hostname)" = "#sk2" ]; then
	sed -i "s/-march=native/-march=skylake/" "Makefile"
fi

config_make_clean
config_make_default_install

config_success
