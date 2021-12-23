#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="numactl"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libnuma.so"
PKG_URL_SRC="numactl-2.0.14.tar.gz"

config_setup

config_package $@

config_configure_make_default_install

config_success
