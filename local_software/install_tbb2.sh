#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="TBB"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libtbb.so"
PKG_URL_SRC="tbb2019_20180718oss_lin.tgz"
PKG_SRC_SUBDIR="tbb2019_20180718oss"

config_setup

config_package $@

echo_info "Installing..."
for i in ./bin ./include ./python ./lib; do
	config_exec cp -r "$i" "$SWEET_LOCAL_SOFTWARE_DST_DIR" || exit 1
done

config_success
