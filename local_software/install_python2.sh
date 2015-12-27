#! /bin/bash

source config.sh

export CC=gcc
export CXX=g++
export LINK=g++


echo "*** Python2 ***"
if [ ! -e "$DST_DIR/bin/python"  -o "$1" != "" ]; then
	SRC_LINK="https://www.python.org/ftp/python/2.7.11/Python-2.7.11.tgz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="Python-2.7.11"

	cd "$SRC_DIR"

	if [ ! -e "$FILENAME" ]; then
		curl "$SRC_LINK" -o "$FILENAME" || exit 1
	fi
	tar xzf "$FILENAME"
	cd "$BASENAME"

	./configure --prefix="$DST_DIR"  || exit 1

	make install || exit 1

	echo "DONE"

else
	echo "Python2 already installed"
fi
