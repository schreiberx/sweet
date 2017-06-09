#! /bin/bash

source config.sh


echo "*** libfreetype ***"

if [ ! -e "$DST_DIR/lib/libfreetype6.so"  -o "$1" != "" ]; then
	SRC_LINK="http://www.martin-schreiber.info/pub/sweet_local_software/freetype-2.8.tar.gz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="freetype-2.8"

	cd "$SRC_DIR"

	if [ ! -e "$FILENAME" ]; then
		curl "$SRC_LINK" -o "$FILENAME" || exit 1
	fi
	tar xzf "$FILENAME"
	cd "$BASENAME"


	./configure --prefix="$DST_DIR" || exit 1
	make install || exit 1

	echo "DONE"

else
	echo "libfreetype2 already installed"
fi
