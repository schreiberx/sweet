#! /bin/bash

source ./install_helpers.sh "" || exit 1
source ./env_vars.sh "" || exit 1

# Name of package
PKG_NAME="gdb"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/gdb"

# URL to source code to fetch it
PKG_URL_SRC="gdb-8.2.tar.xz"

# subdirectory of source in extracted package
# (autodetected with basename of url without file extension if not set)
PKG_SRC_SUBDIR="gdb-8.2"

config_package $@ || exit 1

config_configure_make_default_install || exit 1

config_success || exit 1
