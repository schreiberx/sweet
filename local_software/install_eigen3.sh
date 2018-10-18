#! /bin/bash

source ./install_helpers.sh ""
source ./env_vars.sh ""

PKG_NAME="eigen"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/include/eigen3/Eigen/Eigenvalues"
PKG_URL_SRC="eigen-3.3.3.tar.bz2"

# subdirectory of source in extracted package
# (autodetected with basename of url without file extension if not set)
#PKG_SRC_SUBDIR="eigen-eigen-67e894c6cd8f"

config_setup

config_package $@

mkdir -p build
cd build

echo_info "cmake"
config_exec cmake ../  -DCMAKE_INSTALL_PREFIX="$SWEET_LOCAL_SOFTWARE_DST_DIR" || echo_error_exit "cmake failed"

config_make_default_install

config_success
