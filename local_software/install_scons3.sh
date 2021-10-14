#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="SCons"

if [ "`uname`" == "Darwin" ]; then
	brew install scons
else
	pip3 install scons
fi

config_success

