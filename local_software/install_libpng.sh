#! /bin/bash

source config.sh

echo "*** libpng ***"

if [ ! -e "$DST_DIR/lib/libpng.so"  -o "$1" != "" ]; then
	SRC_LINK="http://www.martin-schreiber.info/pub/sweet_local_software/libpng-1.6.29.tar.gz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="libpng-1.6.29"

	cd "$SRC_DIR"

	if [ ! -e "$FILENAME" ]; then
		curl "$SRC_LINK" -o "$FILENAME" || exit 1
	fi
	tar xzf "$FILENAME"
	cd "$BASENAME"

	# update configure scripts
	./configure --prefix="$DST_DIR" || exit 1
	make install || exit 1

	echo "DONE"

else
	echo "libpng already installed"
fi
