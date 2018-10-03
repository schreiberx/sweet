#! /bin/bash

source ./config.sh ""
source ./env_vars.sh ""


echo "*** FFTW3 ***"
if [ ! -e "$DST_DIR/lib/libfftw3.so"  -o "$1" != "" ]; then
	#SRC_LINK="https://www.martin-schreiber.info/pub/sweet_local_software/fftw-3.3.6-pl2.tar.gz"
	SRC_LINK="https://www.martin-schreiber.info/pub/sweet_local_software/fftw-3.3.8.tar.gz"
	#SRC_LINK="https://www.fftw.org/fftw-3.3.4.tar.gz"
	FILENAME="`basename $SRC_LINK`"
	#BASENAME="fftw-3.3.6-pl2"
	BASENAME="fftw-3.3.8"

	cd "$SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1
	tar xzf "$FILENAME"
	cd "$BASENAME"

	
	CONF_FLAGS=""

	if [ "`uname -s`" == "Linux" -o "`uname -s`" == "Darwin" ]; then
		CONF_FLAGS+=" --enable-openmp "
	fi

	# Activate vectorization code

#	sse only works with single precision
#	CONF_FLAGS+=" --enable-sse"

	CONF_FLAGS+=" --enable-sse2"
	CONF_FLAGS+=" --enable-avx"
	CONF_FLAGS+=" --enable-avx2"
	CONF_FLAGS+=" --enable-avx512"
	CONF_FLAGS+=" --enable-avx-128-fma"

	# Never used directly in Fortran
	CONF_FLAGS+=" --disable-fortran"

	# Enable generation of shared library
	CONF_FLAGS+=" --enable-shared"

	# TODO: Autodetect ARM and use --enable-neon

#	neon only works with single precision
#	CONF_FLAGS+=" --enable-neon"

	echo "Configuration flags: $CONF_FLAGS"

	./configure --prefix="$DST_DIR" $CONF_FLAGS  || exit 1

	make install || exit 1
	make check || exit 1

	echo "DONE"

else
	echo "FFTW already installed"
fi
