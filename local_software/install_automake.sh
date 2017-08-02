#! /bin/bash

source config.sh


#
# Automake
#

echo "*** Automake ***"
if [ "`uname -s`" != "Linux" ]; then
	echo "This script only supports Automake on Linux systems"
else

	if [ ! -e "$DST_DIR/bin/automake"  -o "$1" != "" ]; then
		SRC_LINK="http://www.martin-schreiber.info/pub/sweet_local_software/automake-1.15.tar.gz"
		#SRC_LINK="http://ftp.gnu.org/gnu/automake/automake-1.15.tar.gz"
		FILENAME="`basename $SRC_LINK`"
		BASENAME="automake-1.15"

		cd "$SRC_DIR"

		download "$SRC_LINK" "$FILENAME" || exit 1
		tar xzf "$FILENAME"
		cd "$BASENAME"

		pwd
		./configure --prefix="$DST_DIR" || exit 1
		make -j install

		echo "DONE"

	fi
fi
