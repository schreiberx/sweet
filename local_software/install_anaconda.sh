#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="Anaconda"
PKG_INSTALLED_FILE="$PYTHON_VENV_DIR/bin/activate"
FILENAME="Anaconda3-2020.11-Linux-x86_64.sh"
PKG_URL_SRC="https://repo.anaconda.com/archive/$FILENAME"


config_setup

config_package $@

INSTALLER="$SWEET_LOCAL_SOFTWARE_SRC_DIR/$FILENAME"
sh "${INSTALLER}" -b -u -p ${PYTHON_VENV_DIR} || exit 1

echo
echo "Installation of Anaconda finished"
echo
echo
