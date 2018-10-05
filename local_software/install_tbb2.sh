#! /bin/bash

source ./install_helpers.sh "" || exit 1


# Name of package
PKG_NAME="TBB"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libtbb.so"

# URL to source code to fetch it
PKG_URL_SRC="tbb2019_20180718oss_lin.tgz"

# subdirectory of source in extracted package
# (autodetected with basename of url without file extension if not set)
PKG_SRC_SUBDIR="tbb2019_20180718oss"

config_package $@ || exit 1

echo_info "Installing..."
for i in ./bin ./include ./python ./lib; do
	echo_info " + $i"
	cp -r "$i" "$SWEET_LOCAL_SOFTWARE_DST_DIR" || echo_error_exit "Failed to copy '$i'"
done

config_success || exit 1
