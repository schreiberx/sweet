#! /bin/bash

source config.sh


echo "*** libPFASST ***"
if [ ! -e "$DST_DIR/lib/libpfasst.a"  -o "$1" != "" ]; then
	BASENAME="libpfasst"

	cd "$SRC_DIR"

	if [ ! -e "$BASENAME" ]; then
		git clone https://bitbucket.org/memmett/libpfasst.git || exit 1
	fi

	cd "$BASENAME"

	make || exit 1

	mkdir -p "$DST_DIR/lib/"
	cp -v -f "./lib/libpfasst.a" "$DST_DIR/lib/"
	echo "DONE"

else
	echo "libpfasst (Fortran PFASST) is already installed"
fi
