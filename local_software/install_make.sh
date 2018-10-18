#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="make"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/make"
PKG_URL_SRC="make-4.2.tar.gz"

if [ "`uname -s`" != "Linux" ] && [ "`uname -s`" != "Darwin" ]; then
	echo "This script only supports make on Linux systems"
	exit 1
fi

config_setup

config_package $@

config_configure

config_exec ./build.sh || echo_error_exit "Failed build"

config_exec ./make install || echo_error_exit "Failed to .make install"

config_success
