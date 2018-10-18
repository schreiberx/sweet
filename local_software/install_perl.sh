#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="perl"
PKG_INSTALLED_FILE="${SWEET_LOCAL_SOFTWARE_DST_DIR}/bin/perl"
PKG_URL_SRC="perl-5.28.0.tar.gz"

config_setup

config_package $@

config_exec ./Configure -des -Dprefix=$SWEET_LOCAL_SOFTWARE_DST_DIR -Dnoextensions=ODBM_File

config_make_default_install

config_success
