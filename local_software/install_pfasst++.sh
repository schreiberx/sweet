#! /bin/bash

source ./config.sh ""
source ./env_vars.sh ""


echo "*** PFASST++ ***"
if [ ! -e "$DST_DIR/include/pfasst.hpp"  -o "$1" != "" ]; then
	BASENAME="PFASST"

	cd "$SRC_DIR"

	if [ ! -e "$BASENAME" ]; then
		git clone https://github.com/Parallel-in-Time/PFASST.git || exit 1
	fi

	cd "$BASENAME"

	mkdir -p "$DST_DIR/include/"
	cp -v -f "include/pfasst.hpp" "$DST_DIR/include/"
	cp -v -f -r "include/pfasst" "$DST_DIR/include/"

	echo "DONE"

else
	echo "PFASST++ is already installed"
fi
