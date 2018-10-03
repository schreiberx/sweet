#! /bin/bash

source ./install_helpers.sh "" || exit 1


# Name of package
PKG_NAME="numactl"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libnuma.so"

# URL to source code to fetch it
PKG_URL_SRC="numactl-2.0.12.tar.gz"

config_package $@

config_configure_make_install

config_success
