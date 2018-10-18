#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="likwid"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/likwid-topology"
PKG_URL_SRC="likwid-4.3.2.tar.gz"

config_setup

config_package $@

M_DST="${SWEET_LOCAL_SOFTWARE_DST_DIR}"
M_DST=${M_DST//\//\\/}

sed -i "s/^PREFIX =.*/PREFIX = "${M_DST}"/" config.mk

#sed -i "s/INSTALL_CHOWN = -g root -o root/INSTALL_CHOWN = /" config.mk

sed -i "s/^ACCESSMODE = /&direct#/" config.mk

# Don't build daemon
sed -i "s/^BUILDDAEMON = /&false#/" config.mk

# Don't build setFreq
sed -i "s/^BUILDFREQ = /&false#/" config.mk

config_make_clean

config_make_default_install

config_success
