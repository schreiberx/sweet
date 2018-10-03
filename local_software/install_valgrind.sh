#! /bin/bash

source config.sh ""
source env_vars.sh ""


echo "*** VALGRIND ***"
SRC_LINK="https://www.martin-schreiber.info/pub/sweet_local_software/valgrind-3.13.0.tar.bz2"
FILENAME="`basename $SRC_LINK`"
BASENAME="valgrind-3.13.0"

if [ ! -e "$DST_DIR/bin/valgrind"  -o "$1" != "" ]; then

	cd "$SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1

	tar xjf "$FILENAME"
	cd "$BASENAME"

	
	make clean
	./configure --prefix="$DST_DIR"  || exit 1
	make install || exit 1

	echo "DONE"

else
	echo "SHTNS already installed"
fi
