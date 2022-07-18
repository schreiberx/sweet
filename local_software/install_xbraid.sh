#! /bin/bash

## TODO: rewrite this installation file using functions from install_helpers.sh

source ./install_helpers.sh ""

PKG_NAME="xbraid"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libbraid.a"
PKG_INSTALLED_HEADERS="$SWEET_LOCAL_SOFTWARE_DST_DIR/include/xbraid"
PKG_URL_SRC="https://github.com/XBraid/xbraid.git"

## temporary folder for installation
mkdir -p tmp_xbraid
cd $tmp_xbraid

## download package
git clone $PKG_URL_SRC

## install
cd xbraid
make

## copy library and source files
cp braid/libbraid.a  $PKG_INSTALLED_FILE
mkdir -p $PKG_INSTALLED_HEADERS
cp braid/* $PKG_INSTALLED_HEADERS/.
rm $PKG_INSTALLED_HEADERS/libbraid.a

## clean
cd ../../
rm -r tmp_xbraid




#########config_setup
#########
#########config_package $@
#########
#########
########## Add an empty line to avoid other problems
#########echo "" >> Makefile.local
#########
#########echo_info "Executing 'make clean'..."
#########config_exec make clean
#########
#########config_exec make
#########
#########echo_info "Installing..."
#########
########## Copy modules
#########mkdir -p "$SWEET_LOCAL_SOFTWARE_DST_DIR/include/"
#########echo_info cp -v -f ./include/*mod "$SWEET_LOCAL_SOFTWARE_DST_DIR/include/"
#########cp -v -f ./include/*mod "$SWEET_LOCAL_SOFTWARE_DST_DIR/include/" || echo_error_exit "Failed to install .mod files"
#########
########## Copy static library
#########mkdir -p "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/"
#########echo_info cp -v -f "./lib/libpfasst.a" "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/"
#########cp -v -f "./lib/libpfasst.a" "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/" || echo_error_exit "Failed to install libpfasst.a files"
#########
#########
#########config_success
