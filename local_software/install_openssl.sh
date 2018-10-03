#! /bin/bash

source ./config_install.sh ""
source ./env_vars.sh ""


if [ ! -e "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/libssl.so" -o "$1" != "" ]; then
	SRC_LINK="https://www.martin-schreiber.info/pub/sweet_local_software/openssl-1.1.1.tar.gz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="openssl-1.1.1"

	cd "$SWEET_LOCAL_SOFTWARE_SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1

	tar xzf "$FILENAME"
	cd "$BASENAME"

	./config --prefix="$SWEET_LOCAL_SOFTWARE_DST_DIR" || exit 1
	make
	make install

	echo "DONE"

fi
