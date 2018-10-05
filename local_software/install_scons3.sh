#! /bin/bash

source ./install_helpers.sh "" || exit 1


# Name of package
PKG_NAME="SCons"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/scons"

# URL to source code to fetch it
PKG_URL_SRC="scons-3.0.0.tar.gz"

config_package $@ || exit 1

python3 setup.py install --prefix="$SWEET_LOCAL_SOFTWARE_DST_DIR" || echo_error_exit "Failed to install scons"

config_success || exit 1

