#! /bin/bash

source ./config_install.sh ""
source ./env_vars.sh ""


echo "*** libpng ***"

if [ ! -e "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libpng.so"  -o "$1" != "" ]; then
	SRC_LINK="https://www.martin-schreiber.info/pub/sweet_local_software/libpng-1.6.29.tar.gz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="libpng-1.6.29"

	cd "$SWEET_LOCAL_SOFTWARE_SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1

	tar xzf "$FILENAME"
	cd "$BASENAME"

	# update configure scripts
	./configure --prefix="$SWEET_LOCAL_SOFTWARE_DST_DIR" || exit 1
	make install || exit 1

	echo "DONE"

else
	echo "libpng already installed"
fi
