#! /bin/bash

source ./install_helpers.sh "" || exit 1


# Name of package
PKG_NAME="SHTNS"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libshtns.a"

# URL to source code to fetch it
PKG_URL_SRC="nschaeff-shtns-2018_10_01.tar.bz2"

config_package $@ || exit 1

echo_info_hline
echo_info "SHTNS noOpenMP:"
# Python, no OpenMP
config_configure --disable-openmp || exit 1
config_make_clean || exit 1
config_make_default_install || exit 1

echo_info_hline
echo_info "SHTNS OpenMP:"
# Python, OpenMP
config_configure --enable-openmp || exit 1
config_make_clean || exit 1
config_make_default_install || exit 1

config_success || exit 1
