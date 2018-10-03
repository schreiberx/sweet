#! /bin/bash

source ./config_install.sh ""
source ./env_vars.sh ""


echo "*** RDIC ***"
if [ ! -e "$SWEET_LOCAL_SOFTWARE_DST_DIR/lib/librdic.so.0"  -o "$1" != "" ]; then
	SRC_LINK="https://www.martin-schreiber.info/pub/sweet_local_software/libridc-0.2.tar.gz"
	#SRC_LINK="https://mathgeek.us/files/libridc-0.2.tar.gz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="libridc-0.2"

	cd "$SWEET_LOCAL_SOFTWARE_SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1
	tar xzf "$FILENAME"
	cd "$BASENAME"

	./configure --prefix="$SWEET_LOCAL_SOFTWARE_DST_DIR" || exit 1
	make || exit 1
	make install || exit 1

	echo "DONE"

else
	echo "RDIC is already installed"
fi
