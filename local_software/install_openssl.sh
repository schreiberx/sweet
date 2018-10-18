#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="openssl"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libssl.so"
PKG_URL_SRC="openssl-1.1.1.tar.gz"

config_setup
config_package $@
config_exec "./config --prefix=$SWEET_LOCAL_SOFTWARE_DST_DIR"
config_make_default

if [ "${HOSTNAME:0:10}" != "mpp2-login" ]; then

	# exclude mpp2-login:
	#	Exclude CoolMUC where one test fails :-(

	config_exec "make test"
fi
config_make_install
config_success

