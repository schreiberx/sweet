#! /bin/bash

source config.sh


echo "*** MPACK ***"
if [ -e "$DST_DIR/lib/libmpack.so"  -o "$1" != "" ]; then
	echo "FFTW already installed"
	exit 1
fi


SRC_LINK="https://master.dl.sourceforge.net/project/mplapack/mpack/mpack%200.8.0/mpack-0.8.0-RC2.tar.gz"
FILENAME="`basename $SRC_LINK`"
BASENAME="mpack-0.8.0-RC2.tar.gz"

cd "$SRC_DIR"

curl -C - "$SRC_LINK" -o "$FILENAME" || exit 1

echo "MPACK"
exit 1

tar xzf "$FILENAME"
cd "$BASENAME"

./configure --prefix="$DST_DIR"  || exit 1

echo "DONE"

