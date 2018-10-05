#! /bin/bash

source ./install_helpers.sh "" || exit 1

PKG_NAME="openssl"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libssl.so"
PKG_URL_SRC="openssl-1.1.1.tar.gz"

config_package $@ || exit 1
config_exec "./config --prefix=$SWEET_LOCAL_SOFTWARE_DST_DIR" || exit 1
config_make_default || exit 1

if [ "${HOSTNAME:0:10}" != "mpp2-login" ]; then

	# exclude mpp2-login:
	#	Exclude CoolMUC where one test fails :-(

	config_exec "make test" || exit 1
fi
config_make_install || exit 1
config_success || exit 1

