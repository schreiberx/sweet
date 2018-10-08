#! /bin/bash

source ./install_helpers.sh "" || exit 1

PKG_NAME="git"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/git"
PKG_URL_SRC="git-2.19.0.tar.gz"

if [ "${HOSTNAME:0:10}" != "mpp2-login" ]; then
	echo_error_hline
	echo_error "The git installer is not supported on MPP2"
	echo_error "git pull/push results in empty output without doing anything"
	echo_error_hline
	exit 1
fi

config_package $@ || exit 1
config_configure_make_default || exit 1

config_exec "make test" || exit 1

config_make_install || exit 1

config_success || exit 1

