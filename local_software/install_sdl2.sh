#! /bin/bash

. config.sh


echo "*** SDL2 ***"
if [ "`uname -s`" != "Linux" ]; then
	echo "This script only supports SDL2 on Linux systems"
	echo "Please consider installing pre-build SDL2 packages from http://libsdl.org/release/SDL2-2.0.3.dmg"
else
	if [ ! -e "$DST_DIR/lib/libSDL2.so" ]; then
		SRC_LINK="https://www.libsdl.org/release/SDL2-2.0.3.tar.gz"
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

	fi
fi
