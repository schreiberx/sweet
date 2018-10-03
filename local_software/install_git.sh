#! /bin/bash

source ./config_install.sh ""
source ./env_vars.sh ""


#
# GIT
#


echo "*** GIT ***"
if [ ! -e "$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/git" -o "$1" != "" ]; then
	SRC_LINK="https://www.martin-schreiber.info/pub/sweet_local_software/git-2.19.0.tar.gz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="git-2.19.0"

	cd "$SWEET_LOCAL_SOFTWARE_SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1

	tar xzf "$FILENAME"
	cd "$BASENAME"

	./configure --prefix="$SWEET_LOCAL_SOFTWARE_DST_DIR" --with-openssl || exit 1
	make -j install

	echo "DONE"

fi
