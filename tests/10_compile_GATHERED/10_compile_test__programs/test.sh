#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"

pwd

for i in src/programs/*.cpp; do

	SRC=$(basename $i)
	UNIT_TEST="${SRC/.cpp/}"

	SCONS="scons --gui=disable --sphere-spectral-space=enable --mode=debug"
	SCONS+=" --plane-spectral-space=enable"
	SCONS+=" --sphere-spectral-space=enable"
	SCONS+=" --threading=omp"
	SCONS+=" --eigen=enable"
	SCONS+=" --fortran-source=enable"
	SCONS+=" --sweet-mpi=enable"
	SCONS+=" --libpfasst=enable"
	SCONS+=" --parareal=serial"
	SCONS+=" --program=$UNIT_TEST"

	echo "$SCONS"
	$SCONS  || exit

done
