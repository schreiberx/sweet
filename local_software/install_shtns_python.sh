#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="SHTNS_python"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/python3.6/site-packages/shtns.py"
PKG_URL_SRC="shtns-3.3.1-r694.tar.gz"

config_setup

config_package $@

echo_info_hline
echo_info "SHTNS Python noOpenMP:"
# Python, no OpenMP
config_configure --disable-mem --enable-python --disable-openmp
config_make_clean
config_make_default
python3 setup.py install --prefix="$SWEET_LOCAL_SOFTWARE_DST_DIR" || echo_error_exit "Failed to install"

echo_info_hline
echo_info "SHTNS Python OpenMP:"
# Python, OpenMP
config_configure --disable-mem --enable-python --enable-openmp
config_make_clean
config_make_default
python3 setup.py install --prefix="$SWEET_LOCAL_SOFTWARE_DST_DIR" || echo_error_exit "Failed to install"

config_success
