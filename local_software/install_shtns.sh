#! /bin/bash

source config.sh


echo "*** SHTNS ***"
if [ ! -e "$DST_DIR/lib/libshtns_omp.a"  -o "$1" != "" ]; then
	SRC_LINK="https://bitbucket.org/nschaeff/shtns/downloads/shtns-2.6.6-r535.tar.gz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="shtns-2.6.6-r535"

	cd "$SRC_DIR"

	if [ ! -e "$FILENAME" ]; then
		curl -L "$SRC_LINK" -o "$FILENAME" || exit 1
	fi

	tar xzf "$FILENAME"
	cd "$BASENAME"

	./configure --prefix="$DST_DIR" --enable-openmp || exit 1

	make install || exit 1

	echo "DONE"

else
	echo "SHTNS already installed"
fi
