#! /bin/bash

source ./install_helpers.sh "" || exit 1


# Name of package
PKG_NAME="SDL2"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libSDL2.so"

# URL to source code to fetch it
PKG_URL_SRC="SDL2-2.0.8.tar.gz"

config_package $@

#sed -i -- 's/EXTRA_CFLAGS="$EXTRA_CFLAGS -fpascal-strings"//' ./configure
#./configure --enable-video --enable-video-opengl --prefix="$SWEET_LOCAL_SOFTWARE_DST_DIR" || exit 1

config_configure_make_install

config_success
