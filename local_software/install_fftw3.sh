#! /bin/bash

source config.sh


echo "*** FFTW3 ***"
if [ ! -e "$DST_DIR/lib/libfftw3.so"  -o "$1" != "" ]; then
	SRC_LINK="http://www.martin-schreiber.info/pub/sweet_local_software/fftw-3.3.4.tar.gz"
	#SRC_LINK="http://www.fftw.org/fftw-3.3.4.tar.gz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="fftw-3.3.4"

	cd "$SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1
	tar xzf "$FILENAME"
	cd "$BASENAME"

	
	CONF_FLAGS=""
	if [ "`uname -s`" == "Linux" -o "`uname -s`" == "Darwin" ]; then
		CONF_FLAGS=" --enable-openmp "
	fi

	CONF_FLAGS=" --disable-fortran $CONF_FLAGS"
	./configure --prefix="$DST_DIR" --with-our-malloc16 $CONF_FLAGS --enable-shared  || exit 1
	#./configure $HOST --prefix="$DST_DIR" --with-our-malloc16 $CONF_FLAGS --enable-shared  || exit 1

# AVX is not compiling on all platforms
#	./configure --prefix="$DST_DIR" --with-our-malloc16 $CONF_FLAGS --enable-shared --enable-avx || exit 1

	make install || exit 1

	echo "DONE"

else
	echo "FFTW already installed"
fi
