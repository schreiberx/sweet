#! /bin/bash

source ./install_helpers.sh "" || exit 1

# Name of package
PKG_NAME="cmake"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/cmake"

# URL to source code to fetch it
PKG_URL_SRC="cmake-3.12.3.tar.gz"

config_package $@ || exit 1

config_configure_make_default_install || exit 1

config_make_default_install || exit 1

config_success || exit 1
