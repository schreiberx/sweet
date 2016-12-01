#! /bin/bash

source config.sh


echo "*** PFASST ***"
if [ ! -e "$DST_DIR/lib/librdic.so.0"  -o "$1" != "" ]; then
	BASENAME="PFASST"

	cd "$SRC_DIR"

	if [ ! -e "$BASENAME" ]; then
		git clone https://github.com/Parallel-in-Time/PFASST.git || exit 1
	fi

	cd "$BASENAME"

	cp -v -f "include/pfasst.hpp" "$DST_DIR/include/"
	cp -v -f -r "include/pfasst" "$DST_DIR/include/"

	echo "DONE"

else
	echo "RDIC is already installed"
fi
