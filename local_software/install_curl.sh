#! /bin/bash

source ./install_helpers.sh "" || exit 1


# Name of package
PKG_NAME="curl"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/curl"

# URL to source code to fetch it
PKG_URL_SRC="curl-7.61.1.tar.xz"

config_package $@

config_configure_make_default_install

config_success
