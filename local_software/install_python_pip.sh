#! /bin/bash

source ./install_helpers.sh ""
source ./env_vars.sh ""

PKG_NAME="pip"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/pip"
PKG_URL_SRC="pip-20.0.2.tar.gz"

config_setup

config_package $@

config_exec python3 setup.py build
config_exec python3 setup.py install  --prefix="$SWEET_LOCAL_SOFTWARE_DST_DIR"

config_success
