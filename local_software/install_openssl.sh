#! /bin/bash

source ./install_helpers.sh "" || exit 1

PKG_NAME="openssl"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libssl.so"
PKG_URL_SRC="openssl-1.1.1.tar.gz"

config_package $@
config_exec "./config --prefix=$SWEET_LOCAL_SOFTWARE_DST_DIR"
config_make_default
config_exec "make test"
config_make_install
config_success

