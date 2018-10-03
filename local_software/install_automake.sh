#! /bin/bash

source ./config_install.sh "" || exit 1
source ./env_vars.sh "" || exit 1


# Name of package
PKG_NAME="automake"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/automake"

# URL to source code to fetch it
PKG_URL_SRC="https://www.martin-schreiber.info/pub/sweet_local_software/automake-1.16.tar.xz"

# subdirectory of source in extracted package
# (autodetected with basename of url without file extension if not set)
#SRC_SUBDIR="automake-1.15"

if [ "`uname -s`" != "Linux" ] && [ "`uname -s`" != "Darwin" ]; then
	echo "This script only supports Automake on Linux systems"
	exit 1
fi

config_package $@


echo_info "configure"
./configure --prefix="$SWEET_LOCAL_SOFTWARE_DST_DIR" || echo_error_exit "FAILED configure"

config_make_install

config_success
