#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="valgrind"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/valgrind"
PKG_URL_SRC="valgrind-3.15.0.tar.bz2"

config_setup

config_package $@

config_configure_make_default_install

config_success

