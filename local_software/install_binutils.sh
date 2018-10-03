#! /bin/bash

source ./install_helpers.sh "" || exit 1

echo_error "ERROR: There are some issues if installing binutils such as 'configure: error: memset not found in libc' if configuring libpng\nTherefore, we disable this right here"

# Name of package
PKG_NAME="binutils"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/as"

# URL to source code to fetch it
PKG_URL_SRC="binutils-2.29.tar.gz"

config_package $@

config_configure_make_install

config_success
