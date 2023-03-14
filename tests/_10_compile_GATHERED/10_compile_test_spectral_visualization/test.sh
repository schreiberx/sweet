#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"

if [ "x" != "x$DISPLAY" ]; then 
	echo
	echo "SPECTRAL VISUALIZATION"
	SCONS="scons --program=spectral_visualization --gui=enable --plane-spectral-space=enable --mode=debug"
	echo "$SCONS"
	$SCONS  || exit
fi


mule.benchmark.cleanup_all || exit 1
