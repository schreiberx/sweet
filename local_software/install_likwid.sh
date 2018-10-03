#! /bin/bash

source ./install_helpers.sh "" || exit 1


# Name of package
PKG_NAME="likwid"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/libwid-topology"

# URL to source code to fetch it
PKG_URL_SRC="likwid-4.3.2.tar.gz"

config_package $@

M_SRC="PREFIX = /usr/local"
M_DST="PREFIX = ${SWEET_LOCAL_SOFTWARE_DST_DIR}"
M_SRC=${M_SRC//\//\\/}
M_DST=${M_DST//\//\\/}
sed -i "s/$M_SRC/$M_DST/" config.mk
sed -i "s/INSTALL_CHOWN = -g root -o root/INSTALL_CHOWN = /" config.mk

config_make_install

config_success
