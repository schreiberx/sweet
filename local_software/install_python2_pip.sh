#! /bin/bash

source config.sh

export CC=gcc
export CXX=g++
export LINK=g++

PKG=get-pip.py

echo "*** $PKG ***"
if [ ! -e "$DST_DIR/bin/WHATEVER"  -o "$1" != "" ]; then
	SRC_LINK="https://bootstrap.pypa.io/get-pip.py"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="$PKG"

	cd "$SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1

	python2 get-pip.py || exit 1

	echo "DONE"

else
	echo "$PKG already installed"
fi
