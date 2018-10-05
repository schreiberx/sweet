#! /bin/bash

source ./install_helpers.sh "" || exit 1

# Name of package
PKG_NAME="gcc"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/gcc-8.2"

# URL to source code to fetch it
PKG_URL_SRC="gcc-8.2.0.tar.xz"

config_package $@

echo_info "Downloading additional packages"
./contrib/download_prerequisites || echo_error_exit "Failed to download other required packages"

config_configure --disable-multilib --enable-languages=c++,fortran  --prefix=$SWEET_LOCAL_SOFTWARE_DST_DIR --program-suffix=-8.2

config_make_default_install

for i in g++ gcc gcc-ar gcc-nm gcc-ranlib gfortran gcov gcov-tool gfortran; do
	ln -sf "$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/$i-8.2" "$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/$i" || echo_error_exit "FAILED to create sym links"
done

config_success
