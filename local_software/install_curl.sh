#! /bin/bash

source ./install_helpers.sh "" || exit 1

PKG_NAME="curl"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/curl"
PKG_URL_SRC="curl-7.61.1.tar.xz"

config_package $@ || exit 1

config_configure --with-ca-path=$SWEET_LOCAL_SOFTWARE_DST_DIR/ssl/certs || exit 1

config_make_default || exit 1
config_exec "make test" || exit 1
config_make_install || exit 1
config_success || exit 1

