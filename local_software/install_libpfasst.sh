#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="libpfasst"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libpfasst.a"
PKG_URL_SRC="libpfasst_sweet_2021_11_16.tar.bz2"

config_setup

config_package $@

if [ -z "$MULE_MPIF90" ]; then
	echo_error_exit "MULE_MPIF90 environment variable must be set for libpfasst compilation, see platform configuration"
fi

# Use SWEET's default compiler
sed -i "s/FC.*=.*/FC = ${MULE_MPIF90}/" Makefile.local || echo_error_exit "Replacing Fortran Compiler failed"

# Change flags since there are compiler issues with gfortran 11.x
# We need to do this by inserting a line adding parameters to FFLAGS
sed -i "s/#FFLAGS = -fallow-argument-mismatch/FFLAGS += -fallow-argument-mismatch/" Makefile.local || echo_error_exit "Activating special compiler flag (arg mismatch) failed"

# Deactivate LTO due to bugs in the linker of Ubuntu 21.10
sed -i "1s/#/LDFLAGS = -fno-lto\n#/" Makefile.local || echo_error_exit "Activating special compiler flag (no-lto) failed"

echo_info "Executing 'make clean'..."
config_exec make clean

config_exec make FC=${MULE_MPIF90}

echo_info "Installing..."

# Copy modules
mkdir -p "$SWEET_LOCAL_SOFTWARE_DST_DIR/include/"
echo_info cp -v -f ./include/*mod "$SWEET_LOCAL_SOFTWARE_DST_DIR/include/"
cp -v -f ./include/*mod "$SWEET_LOCAL_SOFTWARE_DST_DIR/include/" || echo_error_exit "Failed to install .mod files"

# Copy static library
mkdir -p "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/"
echo_info cp -v -f "./lib/libpfasst.a" "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/"
cp -v -f "./lib/libpfasst.a" "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/" || echo_error_exit "Failed to install libpfasst.a files"


config_success
