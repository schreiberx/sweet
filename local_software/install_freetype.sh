#! /bin/bash

source config.sh


echo "*** SHTNS ***"
if [ ! -e "$DST_DIR/lib/libfreetype6.so"  -o "$1" != "" ]; then
	SRC_LINK="http://www.martin-schreiber.info/pub/sweet_local_software/freetype-2.7.tar.bz2"
	#SRC_LINK="http://mirror.yannic-bonenberger.com/nongnu/freetype/freetype-2.7.tar.bz2"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="freetype-2.7"

	cd "$SRC_DIR"

	if [ ! -e "$FILENAME" ]; then
		curl "$SRC_LINK" -o "$FILENAME" || exit 1
	fi

	tar xjf "$FILENAME"
	cd "$BASENAME"

	# no OpenMP
	make clean
	./configure --prefix="$DST_DIR" || exit 1
	make install || exit 1

	echo "DONE"

else
	echo "LIBFREETYPE already installed"
fi
