#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="libfreetype"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libfreetype.so"
PKG_URL_SRC="freetype-2.8.tar.gz"

config_setup

config_package $@

config_configure_make_default_install

config_success
