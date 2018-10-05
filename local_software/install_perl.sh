#! /bin/bash

source ./install_helpers.sh "" || exit 1

# Name of package
PKG_NAME="perl"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="${SWEET_LOCAL_SOFTWARE_DST_DIR}/bin/perl"

# URL to source code to fetch it
PKG_URL_SRC="perl-5.28.0.tar.gz"

config_package $@

config_exec ./Configure -des -Dprefix=$SWEET_LOCAL_SOFTWARE_DST_DIR -Dnoextensions=ODBM_File

config_make_default_install

config_success
