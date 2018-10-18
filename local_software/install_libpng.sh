#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="libpng"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libpng.so"
PKG_URL_SRC="libpng-1.6.29.tar.gz"

config_setup

config_package $@

config_configure_make_default_install

config_success
