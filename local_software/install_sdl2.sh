#! /bin/bash

source config.sh


echo "*** SDL2 ***"

if [ ! -e "$DST_DIR/lib/libSDL2.so"  -o "$1" != "" ]; then
	SRC_LINK="http://www.martin-schreiber.info/pub/sweet_local_software/SDL2-2.0.3.tar.gz"
	#SRC_LINK="https://www.libsdl.org/release/SDL2-2.0.3.tar.gz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="SDL2-2.0.3"

	cd "$SRC_DIR"

	if [ ! -e "$FILENAME" ]; then
		curl "$SRC_LINK" -o "$FILENAME" || exit 1
	fi
	tar xzf "$FILENAME"
	cd "$BASENAME"

	# update configure scripts
	#sh autogen.sh
	sed -i -- 's/EXTRA_CFLAGS="$EXTRA_CFLAGS -fpascal-strings"//' ./configure
	./configure --prefix="$DST_DIR" || exit 1
	make install || exit 1

	echo "DONE"

else
	echo "SDL2 already installed"
fi
