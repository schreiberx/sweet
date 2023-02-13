#! /bin/bash

#
# OpenMPI is deprecated since this is troublesome with valgrind
#

source ./install_helpers.sh ""

PKG_NAME="OpenMPI"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/openmpi.so"
PKG_URL_SRC="openmpi-4.1.4.tar.bz2"

config_setup

config_package $@

config_configure --enable-mpi-fortran

config_make_default_install

config_success
