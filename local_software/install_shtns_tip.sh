#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="SHTNS"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libshtns.a"
PKG_URL_SRC="https://bitbucket.org/nschaeff/shtns/get/tip.tar.bz2"

config_setup

config_package $@

echo_info_hline
echo_info "SHTNS noOpenMP:"
config_configure --disable-mem --disable-openmp
#config_configure --disable-openmp
config_make_clean
config_make_default_install

echo_info_hline
echo_info "SHTNS OpenMP:"
config_configure --disable-mem --enable-openmp
#config_configure --enable-openmp
config_make_clean
config_make_default_install

config_success
