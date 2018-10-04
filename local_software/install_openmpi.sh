#! /bin/bash

source ./install_helpers.sh "" || exit 1


# Name of package
PKG_NAME="OpenMPI"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/mpicc"

# URL to source code to fetch it
PKG_URL_SRC="openmpi-3.1.2.tar.bz2"

config_package $@

config_configure --enable-mpi-fortran
config_make_default_install

config_success
