#! /bin/bash

source config.sh


echo "*** BINUTILS ***"
if [ -e "$DST_DIR/bin/as" -a "$1" == "" ]; then
	echo "Binutils already installed"
	exit 1
fi


SRC_LINK="https://www.martin-schreiber.info/pub/sweet_local_software/binutils-2.29.tar.gz"
FILENAME="`basename $SRC_LINK`"
BASENAME="binutils-2.29"

cd "$SRC_DIR" || exit 1

download "$SRC_LINK" "$FILENAME" || exit 1
tar xzf "$FILENAME" || exit 1
cd "$BASENAME" || exit 1

./configure --prefix="$DST_DIR" || exit 1
make || exit 1
make install || exit 1

echo "DONE"
