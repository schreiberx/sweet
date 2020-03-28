#! /bin/bash

source ./install_helpers.sh ""

# Setup environment
#source "$MULE_SOFTWARE_ROOT/local_software/local/python_env/bin/activate"

PKG_NAME="pip packages"

#echo ln -sf "$MULE_SOFTWARE_ROOT/local_software/local/python_env/bin/pip" "$MULE_SOFTWARE_ROOT/local_software/local/python_env/bin/pip3"

echo_info "Upgrading pip..."
pip3 install --upgrade --user pip -v

PIP_PACKAGES="numpy matplotlib scipy"
echo_info "Installing $PIP_PACKAGES..."
pip3 install --user $PIP_PACKAGES -v || echo_error_exit "Failed to install $PIP_PACKAGES"

config_success
