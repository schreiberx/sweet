#! /bin/bash

source config.sh

#export CC=gcc
#export CXX=g++
#export LINK=g++


echo "*** Python3 ***"

if [ ! -e "$DST_DIR/bin/python3"  -o "$1" != "" ]; then
	SRC_LINK="https://www.martin-schreiber.info/pub/sweet_local_software/Python-3.6.2.tgz"
	#SRC_LINK="https://www.python.org/ftp/python/2.7.11/Python-2.7.11.tgz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="Python-3.6.2"

	cd "$SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1
	tar xzf "$FILENAME"
	cd "$BASENAME"

	#sed -i "s/#zlib/zlib/" Modules/Setup	 || exit 1
	#sed -i "s/#binascii/binascii/" Modules/Setup	|| exit 1

	./configure --prefix="$DST_DIR" || exit 1

	make install || exit 1
	ln -sf "$DST_DIR/bin/python3" "$DST_DIR/bin/python" 

	echo "DONE"

else
	echo "Python3 already installed"
fi
