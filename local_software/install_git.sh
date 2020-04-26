#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="git"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/git"
PKG_URL_SRC="git-2.19.0.tar.gz"

if [ "${HOSTNAME:0:10}" != "mpp2-login" ]; then
	echo_error_hline
	echo_error "The git installer is only supported on MPP2"
	echo_error "git pull/push results in empty output without doing anything"
	echo_error_hline
	exit 1
fi

config_setup

config_package $@
config_configure_make_default

config_exec make test

config_make_install

config_success

