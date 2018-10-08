#! /bin/bash

source ./install_helpers.sh "" || exit 1


# Name of package
PKG_NAME="Python3"

# Path to one file of installed package to test for existing installation
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/python3"

# URL to source code to fetch it
PKG_URL_SRC="Python-3.6.2.tgz"

config_package $@ || exit 1

config_configure_make_default_install || exit 1

#mkdir -p "$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/python"
#Create symlink to make python linked to python3
config_exec ln -sf "$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/python3" "$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/python" 

echo_info "Upgrading pip..."
pip3 install --upgrade pip

echo_info "Installing numpy and matplotlib..."
pip3 install numpy matplotlib || echo_error_exit "Failed to install numpy and matplotlib"

config_success || exit 1
