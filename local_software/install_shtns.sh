#! /bin/bash

source ./install_helpers.sh "" || exit 1


# Name of package
PKG_NAME="SHTNS Python"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libshtns.a"

# URL to source code to fetch it
PKG_URL_SRC="nschaeff-shtns-2018_10_01.tar.bz2"

config_package $@

echo_info_hline
echo_info "SHTNS noOpenMP:"
# Python, no OpenMP
config_configure --disable-openmp
config_make_clean
config_make_install

echo_info_hline
echo_info "SHTNS OpenMP:"
# Python, OpenMP
make_configure --enable-openmp
config_make_clean
config_make_install

config_success
