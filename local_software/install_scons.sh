#! /bin/bash

source config.sh


echo "*** SCONS ***"
if [ ! -e "$DST_DIR/bin/scons"  -o "$1" != "" ]; then
	SRC_LINK="http://kent.dl.sourceforge.net/project/scons/scons/2.4.0/scons-2.4.0.tar.gz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="scons-2.4.0"

	cd "$SRC_DIR"

	if [ ! -e "$FILENAME" ]; then
		curl "$SRC_LINK" -o "$FILENAME" || exit 1
	fi
	tar xzf "$FILENAME"
	cd "$BASENAME"

	python setup.py install --prefix="$DST_DIR" || exit 1

	echo "DONE"

fi
