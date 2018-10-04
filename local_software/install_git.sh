#! /bin/bash

source ./install_helpers.sh "" || exit 1


# Name of package
PKG_NAME="git"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/git"

# URL to source code to fetch it
PKG_URL_SRC="git-2.19.0.tar.gz"

config_package $@

config_configure_make_default_install

config_success
