#! /bin/bash

source config.sh


echo "*** SCONS ***"
if [ ! -e "$DST_DIR/bin/scons"  -o "$1" != "" ]; then
	SRC_LINK="http://www.martin-schreiber.info/pub/sweet_local_software/scons-2.4.0.tar.gz"
	#SRC_LINK="https://netix.dl.sourceforge.net/project/scons/scons/2.4.0/scons-2.4.0.tar.gz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="scons-2.4.0"

	cd "$SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1
	tar xzf "$FILENAME" || exit 1
	cd "$BASENAME" || exit 1

	python setup.py install --prefix="$DST_DIR" || exit 1

	echo "DONE"

else
	echo "SCONS is already installed"
fi
