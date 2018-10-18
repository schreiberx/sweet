#! /bin/bash

cd "$SWEET_ROOT"

if [ "x" != "x$DISPLAY" ]; then 
	echo
	echo "SPECTRAL VISUALIZATION"
	SCONS="scons --program=spectral_visualization --gui=enable --plane-spectral-space=enable --mode=debug"
	echo "$SCONS"
	$SCONS  || exit
fi


