#! /bin/bash

source config.sh

export CC=gcc
export CXX=g++
export LINK=g++


echo "*** mpmath ***"
if [ ! -e "$DST_DIR/bin/WHATEVER"  -o "$1" != "" ]; then
	SRC_LINK="http://www.martin-schreiber.info/pub/sweet_local_software/mpmath-0.19.tar.gz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="mpmath-0.19"

	cd "$SRC_DIR"

	if [ ! -e "$FILENAME" ]; then
		curl "$SRC_LINK" -o "$FILENAME" || exit 1
	fi
	tar xzf "$FILENAME"
	cd "$BASENAME"

	python2 setup.py build || exit 1
	python2 setup.py install --prefix="$DST_DIR" || exit 1

	echo "DONE"

else
	echo "mpmath already installed"
fi
