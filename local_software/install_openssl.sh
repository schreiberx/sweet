#! /bin/bash

source ./install_helpers.sh "" || exit 1


# Name of package
PKG_NAME="openssl"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libssl.so"

# URL to source code to fetch it
PKG_URL_SRC="openssl-1.1.1.tar.gz"

config_package $@

M="./config --prefix=$SWEET_LOCAL_SOFTWARE_DST_DIR $@"
echo_info "Executing '${M}'"
$M || echo_error_exit "FAILED '${M}'"

config_make_default_install

config_success

