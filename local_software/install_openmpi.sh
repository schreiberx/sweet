#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="OpenMPI"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/mpicc"
PKG_URL_SRC="openmpi-3.1.2.tar.bz2"

config_setup

config_package $@

config_configure --enable-mpi-fortran

config_make_default_install

config_success
