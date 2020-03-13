#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="Python3"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/python3"
PKG_URL_SRC="Python-3.8.2.tar.xz"

config_setup

config_package $@

config_configure_make_default_install

#mkdir -p "$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/python"
#Create symlink to make python linked to python3
#config_exec ln -sf "$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/python3" "$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/python" 

echo_info "Upgrading pip..."
pip3 install --upgrade pip

PIP_PACKAGES="numpy matplotlib scipy"
echo_info "Installing $PIP_PACKAGES..."
pip3 install $PIP_PACKAGES || echo_error_exit "Failed to install $PIP_PACKAGES"

config_success
