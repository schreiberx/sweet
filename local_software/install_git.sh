#! /bin/bash

source config.sh
source env_vars.sh


#
# GIT
#


echo "*** GIT ***"
if [ ! -e "$DST_DIR/bin/git" -o "$1" != "" ]; then
	SRC_LINK="https://www.martin-schreiber.info/pub/sweet_local_software/git-2.19.0.tar.gz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="git-2.19.0"

	cd "$SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1

	tar xzf "$FILENAME"
	cd "$BASENAME"

	./configure --prefix="$DST_DIR" || exit 1
	make -j install

	echo "DONE"

fi
