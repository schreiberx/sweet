#! /bin/bash

source ./install_helpers.sh ""

PKG_NAME="SCons"


if [ "${MULE_PLATFORM_ID:0:9}" == "supermuc_" ]; then
	echo_error ""
	echo_error "Can't execute:"
	echo_error "	pip3 install scons==4.0"
	echo_error ""
	echo_error "You need to install 'scons' manually since there's no direct internet connection!"
	echo_error ""
	exit 1
fi

if [ "`uname`" == "Darwin" ]; then
	brew install scons
else
	pip3 install scons
fi

config_success

