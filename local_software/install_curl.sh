#! /bin/bash

source ./install_helpers.sh "" || exit 1

PKG_NAME="curl"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/curl"
PKG_URL_SRC="curl-7.61.1.tar.xz"

config_package $@
config_configure_make_default
config_exec "make test"
config_configure_make_install
config_success

