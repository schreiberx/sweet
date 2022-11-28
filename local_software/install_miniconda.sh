#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="Miniconda"
PKG_INSTALLED_FILE="$PYTHON_VENV_DIR/bin/activate"
FILENAME="Miniconda3-py39_4.12.0-Linux-x86_64.sh"
PKG_URL_SRC="https://www.martin-schreiber.info/pub/sweet_local_software/$FILENAME"

if [ "`uname`" == "Darwin" ]; then
	echo_error_hline
	echo_error "Miniconda is not required on MacOSX with homebrew, see INSTALL_MACOSX"
	echo_error_hline
	return
fi


config_setup

config_package $@

INSTALLER="$SWEET_LOCAL_SOFTWARE_SRC_DIR/$FILENAME"
sh "${INSTALLER}" -b -u -p ${PYTHON_VENV_DIR} || exit 1

. ${PYTHON_VENV_DIR}/bin/activate || exit 1

pip3 install matplotlib numpy sympy scipy

echo
echo "Installation of Miniconda finished"
echo
echo
