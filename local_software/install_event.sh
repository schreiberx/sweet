#! /bin/bash

source ./install_helpers.sh ""


PKG_NAME="libevent"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libevent.a"
PKG_URL_SRC="libevent-2.1.11-stable.tar.gz"

config_setup

config_package $@

config_configure_make_default_install

config_success


