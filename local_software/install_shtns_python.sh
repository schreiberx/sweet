#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="SHTNS_python"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/python3.6/site-packages/shtns.py"

PKG_URL_SRC="shtns-3.3.1-r694.tar.gz"

config_setup

config_package $@

if [ "#$TRAVIS" != "#" ]; then
	echo_info "Detected Travis"
	echo_info "Disabling SIMD because of Travis"
	EXTRA_FLAGS="--disable-mkl --disable-knl --disable-cuda --disable-simd"
fi


echo_info_hline
echo_info "SHTNS Python noOpenMP:"
# Python, no OpenMP
config_configure --enable-python --disable-openmp $EXTRA_FLAGS
config_make_clean
config_make_default
python3 setup.py install --prefix="$SWEET_LOCAL_SOFTWARE_DST_DIR" || echo_error_exit "Failed to install"

echo_info_hline
echo_info "SHTNS Python OpenMP:"
# Python, OpenMP
config_configure --enable-python --enable-openmp $EXTRA_FLAGS
config_make_clean
config_make_default
python3 setup.py install --prefix="$SWEET_LOCAL_SOFTWARE_DST_DIR" || echo_error_exit "Failed to install"

config_success
