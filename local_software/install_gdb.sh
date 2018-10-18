#! /bin/bash

source ./install_helpers.sh ""
source ./env_vars.sh ""

PKG_NAME="gdb"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/gdb"
PKG_URL_SRC="gdb-8.2.tar.xz"

config_setup

config_package $@

config_configure_make_default_install

config_success
