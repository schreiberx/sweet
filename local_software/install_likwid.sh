#! /bin/bash

source ./config_install.sh ""
source ./env_vars.sh ""


echo "*** LIKWID ***"
if [ ! -e "$SWEET_LOCAL_SOFTWARE_DST_DIR/bin/libwid-topology" -o "$1" != "" ]; then
	SRC_LINK="https://ftp.fau.de/pub/likwid/likwid-4.2.0.tar.gz"
	FILENAME="`basename $SRC_LINK`"
	BASENAME="likwid-4.2.0"

	cd "$SWEET_LOCAL_SOFTWARE_SRC_DIR"

	download "$SRC_LINK" "$FILENAME" || exit 1

	tar xzf "$FILENAME"
	cd "$BASENAME"

	M_SRC="PREFIX = /usr/local"
	M_DST="PREFIX = ${SWEET_LOCAL_SOFTWARE_DST_DIR}"
	M_SRC=${M_SRC//\//\\/}
	M_DST=${M_DST//\//\\/}
	sed -i "s/$M_SRC/$M_DST/" config.mk
	sed -i "s/INSTALL_CHOWN = -g root -o root/INSTALL_CHOWN = /" config.mk

	#./configure --prefix="$SWEET_LOCAL_SOFTWARE_DST_DIR" || exit 1
	make || exit 1
	make install || exit 1

	echo "DONE"

else
	echo "LIKWID already installed"
fi
