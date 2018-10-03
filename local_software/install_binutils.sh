#! /bin/bash

source ./config_install.sh "" || exit 1
source ./env_vars.sh "" || exit 1


# Name of package
PKG_NAME="binutils"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/as"

# URL to source code to fetch it
PKG_URL_SRC="https://www.martin-schreiber.info/pub/sweet_local_software/binutils-2.29.tar.gz"

# subdirectory of source in extracted package
# (autodetected with basename of url without file extension if not set)
#SRC_SUBDIR="automake-1.15"

config_package $@

echo_info "configure"
./configure --prefix="$SWEET_LOCAL_SOFTWARE_DST_DIR" || echo_error_exit "FAILED configure"

config_make_install

config_success
