#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="SDL2"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libSDL2.so"
PKG_URL_SRC="SDL2-2.0.8.tar.gz"

config_setup

config_package $@

config_configure_make_default_install

config_success
