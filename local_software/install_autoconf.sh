#! /bin/bash

source config.sh


#
# Autoconf
#

echo "*** Autoconf ***"
if [ "`uname -s`" != "Linux" ] && [ "`uname -s`" != "Darwin" ]; then
	echo "This script only supports Autoconf on Linux systems"
else

	if [ ! -e "$DST_DIR/bin/autoconf" -o "$1" != ""  ]; then
		SRC_LINK="http://www.martin-schreiber.info/pub/sweet_local_software/autoconf-2.69.tar.gz"
		#SRC_LINK="http://ftp.gnu.org/gnu/autoconf/autoconf-2.69.tar.gz"
		FILENAME="`basename $SRC_LINK`"
		BASENAME="autoconf-2.69"

		cd "$SRC_DIR"

		if [ ! -e "$FILENAME" ]; then
			curl "$SRC_LINK" -o "$FILENAME" || exit 1
		fi
		tar xzf "$FILENAME"
		cd "$BASENAME"

		pwd
		./configure --prefix="$DST_DIR" || exit 1
		make -j install

		echo "DONE"

	else
		echo "Autoconf already installed"
	fi
fi
