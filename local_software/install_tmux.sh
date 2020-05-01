#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="tmux"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/tmux"
PKG_URL_SRC="tmux-3.1a.tar.gz"

config_setup

config_package $@

config_configure_make_default_install

config_success
