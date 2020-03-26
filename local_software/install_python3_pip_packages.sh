#! /bin/bash

source ./install_helpers.sh ""

echo_info "Upgrading pip..."
pip3 install --upgrade pip

PIP_PACKAGES="numpy matplotlib scipy"
echo_info "Installing $PIP_PACKAGES..."
pip3 install $PIP_PACKAGES || echo_error_exit "Failed to install $PIP_PACKAGES"

config_success
