#! /bin/bash

source ./install_helpers.sh "" || exit 1

PKG_NAME="git"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/git"
PKG_URL_SRC="git-2.19.0.tar.gz"

config_package $@ || exit 1
config_configure_make_default || exit 1

if [ "${HOSTNAME:0:10}" != "mpp2-login" ]; then
	# exclude mpp2-login:
	#	tests successful, but removal of test directory fails

	config_exec "make test" || exit 1
fi

config_make_install || exit 1

config_success || exit 1

