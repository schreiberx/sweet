#! /bin/bash

source ./config_install.sh "" || exit 1
source ./env_vars.sh "" || exit 1


# Name of package
PKG_NAME="autoconf"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/autoconf"

# URL to source code to fetch it
PKG_URL_SRC="https://www.martin-schreiber.info/pub/sweet_local_software/autoconf-2.69.tar.gz"

# subdirectory of source in extracted package
# (autodetected with basename of url without file extension if not set)
#SRC_SUBDIR="autoconf-2.69"

if [ "`uname -s`" != "Linux" ] && [ "`uname -s`" != "Darwin" ]; then
	echo "This script only supports Autoconf on Linux systems"
	exit 1
fi

config_package $@

#
# We are now in the extracted source folder
#

echo_info "configure"
./configure --prefix="$SWEET_LOCAL_SOFTWARE_DST_DIR" || echo_error_exit "FAILED configure"

config_make_install

config_success
