#! /bin/bash

source config.sh


echo "*** RDIC ***"
if [ ! -e "$DST_DIR/lib/librdic.so.0"  -o "$1" != "" ]; then
	SRC_LINK="http://mathgeek.us/files/libridc-0.2.tar.gz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="libridc-0.2"

	cd "$SRC_DIR"

	if [ ! -e "$FILENAME" ]; then
		curl "$SRC_LINK" -o "$FILENAME" || exit 1
	fi
	tar xzf "$FILENAME"
	cd "$BASENAME"

	./configure --prefix="$DST_DIR" || exit 1
	make || exit 1
	make install || exit 1

	echo "DONE"

else
	echo "RDIC is already installed"
fi
