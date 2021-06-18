#! /bin/bash

source ./install_helpers.sh ""

# Override this for this install script since this is troublesome on the LRZ Linuxcluster
TMPDIR="$SWEET_LOCAL_SOFTWARE_SRC_DIR/tmp/"
mkdir -p "$TMPDIR"

PKG_NAME="Python3"
PKG_INSTALLED_FILE="$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/python3"
PKG_URL_SRC="Python-3.8.2.tar.xz"


# Remove local python environment directory
rm -f "$MULE_SOFTWARE_ROOT/local_software/local/python_env"


config_setup

config_package $@

config_configure_make_default_install

#mkdir -p "$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/python"
#Create symlink to make python linked to python3
#config_exec ln -sf "$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/python3" "$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/python" 

# Use SWEET python environment in case that the system-wide installed python is used
# Ignore errors in case that this folder doesn't exist
python3 -m venv "$SCRIPTDIR/local/python_env"

# Setup environment
source "$MULE_SOFTWARE_ROOT/local_software/local/python_env/bin/activate"


config_success
