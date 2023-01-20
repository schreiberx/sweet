#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="git"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/git"
PKG_URL_SRC="git-2.39.0.tar.gz"

config_setup

config_package $@

# Use --with-curl to support https
config_configure_make_default --with-curl

config_exec make test

config_make_install

config_success

