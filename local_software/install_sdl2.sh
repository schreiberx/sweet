#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="SDL2"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libSDL2.so"
PKG_URL_SRC="SDL2-2.26.2.tar.gz"

config_setup

config_package $@

config_configure --enable-video-opengl
config_make_default_install

config_success
