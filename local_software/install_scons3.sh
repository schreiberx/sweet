#! /bin/bash

source config.sh
source env_vars.sh


echo "*** SCONS ***"
if [ ! -e "$DST_DIR/bin/scons"  -o "$1" != "" ]; then
	SRC_LINK="https://www.martin-schreiber.info/pub/sweet_local_software/scons-3.0.0.tar.gz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="scons-3.0.0"

	cd "$SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1
	tar xzf "$FILENAME" || exit 1
	cd "$BASENAME" || exit 1

	python3 setup.py install --prefix="$DST_DIR" || exit 1

	echo "DONE"

else
	echo "SCONS is already installed"
fi
