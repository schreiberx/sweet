#! /bin/bash

source ./config_install.sh ""
source ./env_vars.sh ""


echo "*** PFASST++ ***"
if [ ! -e "$SWEET_LOCAL_SOFTWARE_DST_DIR/include/pfasst.hpp"  -o "$1" != "" ]; then
	BASENAME="PFASST"

	cd "$SWEET_LOCAL_SOFTWARE_SRC_DIR"

	if [ ! -e "$BASENAME" ]; then
		git clone https://github.com/Parallel-in-Time/PFASST.git || exit 1
	fi

	cd "$BASENAME"

	mkdir -p "$SWEET_LOCAL_SOFTWARE_DST_DIR/include/"
	cp -v -f "include/pfasst.hpp" "$SWEET_LOCAL_SOFTWARE_DST_DIR/include/"
	cp -v -f -r "include/pfasst" "$SWEET_LOCAL_SOFTWARE_DST_DIR/include/"

	echo "DONE"

else
	echo "PFASST++ is already installed"
fi
