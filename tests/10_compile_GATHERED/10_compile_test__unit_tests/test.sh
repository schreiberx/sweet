#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"

pwd

for i in src/unit_tests/*.cpp; do

	SRC=$(basename $i)
	UNIT_TEST="${SRC/.cpp/}"

	SCONS="scons --gui=disable --sphere-spectral-space=enable --mode=debug"
	SCONS+=" --plane-spectral-space=enable"
	SCONS+=" --sphere-spectral-space=enable"
	SCONS+=" --threading=omp"
	SCONS+=" --eigen=enable"
	SCONS+=" --parareal=serial"
	SCONS+=" --parareal-sphere=enable"
	SCONS+=" --unit-test=$UNIT_TEST"

	echo "$SCONS"
	$SCONS  || exit

done
