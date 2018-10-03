#! /bin/bash

source ./config_install.sh ""
source ./env_vars.sh ""


echo "*** libfreetype ***"

if [ ! -e "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libfreetype6.so"  -o "$1" != "" ]; then
	SRC_LINK="https://www.martin-schreiber.info/pub/sweet_local_software/freetype-2.8.tar.gz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="freetype-2.8"

	cd "$SWEET_LOCAL_SOFTWARE_SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1

	tar xzf "$FILENAME"
	cd "$BASENAME"


	./configure --prefix="$SWEET_LOCAL_SOFTWARE_DST_DIR" || exit 1
	make install || exit 1

	echo "DONE"

else
	echo "libfreetype2 already installed"
fi
