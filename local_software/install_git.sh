#! /bin/bash

source ./install_helpers.sh "" || exit 1

PKG_NAME="git"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/git"
PKG_URL_SRC="git-2.19.0.tar.gz"

config_package $@
config_configure_make_default
config_exec "make test"
config_make_install
config_success

