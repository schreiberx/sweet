#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="SCons"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/scons"
PKG_URL_SRC="scons-3.1.2.tar.gz"

config_setup

config_package $@

config_exec python3 setup.py install --prefix="$SWEET_LOCAL_SOFTWARE_DST_DIR" || echo_error_exit "Failed to install scons"

sed -i "s/env python/env python3/" "$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/scons"

config_success

