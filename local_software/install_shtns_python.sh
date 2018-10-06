#! /bin/bash

source ./install_helpers.sh "" || exit 1


# Name of package
PKG_NAME="SHTNS Python"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/python3.6/site-packages/shtns.py"

# URL to source code to fetch it
PKG_URL_SRC="nschaeff-shtns-2018_10_01.tar.bz2"

config_package $@ || exit 1

echo_info_hline
echo_info "SHTNS Python noOpenMP:"
# Python, no OpenMP
config_configure --enable-python --disable-openmp || exit 1
config_make_clean || exit 1
config_make_default || exit 1
python3 setup.py install --prefix="$SWEET_LOCAL_SOFTWARE_DST_DIR" || echo_error_exit "Failed to install"

echo_info_hline
echo_info "SHTNS Python OpenMP:"
# Python, OpenMP
config_configure --enable-python --enable-openmp || exit 1
config_make_clean || exit 1
config_make_default || exit 1
python3 setup.py install --prefix="$SWEET_LOCAL_SOFTWARE_DST_DIR" || echo_error_exit "Failed to install"

config_success || exit 1
