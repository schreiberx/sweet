#! /bin/bash

source ./install_helpers.sh "" || exit 1
source ./env_vars.sh "" || exit 1

# Name of package
PKG_NAME="eigen"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/include/eigen3/Eigen/Eigenvalues"

# URL to source code to fetch it
PKG_URL_SRC="eigen-3.3.3.tar.bz2"

# subdirectory of source in extracted package
# (autodetected with basename of url without file extension if not set)
PKG_SRC_SUBDIR="eigen-eigen-67e894c6cd8f"

config_package $@

mkdir -p build
cd build

echo_info "cmake"
cmake ../  -DCMAKE_INSTALL_PREFIX="$SWEET_LOCAL_SOFTWARE_DST_DIR" || echo_error_exit "cmake failed"

config_make_install

config_success
