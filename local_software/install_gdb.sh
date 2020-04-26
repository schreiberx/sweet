#! /bin/bash

source ./install_helpers.sh ""
source ./env_vars.sh ""

PKG_NAME="gdb"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/gdb"
PKG_URL_SRC="gdb-9.1.tar.xz"

config_setup

config_package $@

config_exec mkdir -p build
config_exec cd build
config_exec ../configure --prefix=$SWEET_LOCAL_SOFTWARE_DST_DIR
config_make_default
config_make_install
#config_configure_make_default_install

config_success
