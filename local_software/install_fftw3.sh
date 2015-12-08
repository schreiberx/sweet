#! /bin/bash

source config.sh


echo "*** FFTW3 ***"
if [ ! -e "$DST_DIR/lib/libfftw3.so" ]; then
	SRC_LINK="http://www.fftw.org/fftw-3.3.4.tar.gz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="fftw-3.3.4"

	cd "$SRC_DIR"

	if [ ! -e "$FILENAME" ]; then
		curl "$SRC_LINK" -o "$FILENAME" || exit 1
	fi
	tar xzf "$FILENAME"
	cd "$BASENAME"

	
	CONF_FLAGS=""
	if [ "`uname -s`" == "Linux" ]; then
		CONF_FLAGS=" --enable-openmp "
	fi
	./configure --prefix="$DST_DIR" --with-our-malloc16 $CONF_FLAGS --enable-shared  || exit 1

# AVX is not compiling on all platforms
#	./configure --prefix="$DST_DIR" --with-our-malloc16 $CONF_FLAGS --enable-shared --enable-avx || exit 1

	make install || exit 1

	echo "DONE"

fi
