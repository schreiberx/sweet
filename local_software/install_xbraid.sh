#! /bin/bash


source ./install_helpers.sh ""

PKG_NAME="xbraid"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libbraid.a"
PKG_INSTALLED_HEADERS="$SWEET_LOCAL_SOFTWARE_DST_DIR/include/xbraid"
PKG_URL_SRC="https://github.com/XBraid/xbraid/archive/refs/tags/v3.1.0.tar.gz"

config_setup

config_package $@

# Add an empty line to avoid other problems
echo "" >> Makefile.local

echo_info "Executing 'make clean'..."
config_exec make clean

config_exec make braid

echo_info "Installing..."

# Copy modules
mkdir -p $PKG_INSTALLED_HEADERS
echo_info cp -v -f braid/* $PKG_INSTALLED_HEADERS/.
cp -v -f  braid/* $PKG_INSTALLED_HEADERS/. || echo_error_exit "Failed to install .mod files"

# Copy static library
mkdir -p "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/"
echo_info cp -v -f braid/libbraid.a $PKG_INSTALLED_FILE
cp -v -f braid/libbraid.a $PKG_INSTALLED_FILE || echo_error_exit "Failed to install libbraid.a files"

config_success
