#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="cmake"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/cmake"
PKG_URL_SRC="cmake-3.12.3.tar.gz"

config_setup

config_package $@

config_configure_make_default_install

config_make_default_install

config_success
