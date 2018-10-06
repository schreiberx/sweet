#! /bin/bash

source ./install_helpers.sh "" || exit 1

# Name of package
PKG_NAME="make"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/make"

# URL to source code to fetch it
PKG_URL_SRC="make-4.2.tar.gz"

if [ "`uname -s`" != "Linux" ] && [ "`uname -s`" != "Darwin" ]; then
	echo "This script only supports make on Linux systems"
	exit 1
fi

config_package $@ || exit 1

config_configure || exit 1

./build.sh || echo_error_exit "Failed build"

./make install || echo_error_exit "Failed to .make install"

config_success || exit 1
