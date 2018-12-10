#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"

for i in 0 1 2 3; do
	SCONS="scons --program=advection --numa-block-allocator=$i --threading=omp --mode=debug"
	echo "$SCONS"
	$SCONS || exit

	SCONS="scons --program=advection --numa-block-allocator=$i --threading=off --mode=debug"
	echo "$SCONS"
	$SCONS || exit

	SCONS="scons --program=advection --numa-block-allocator=$i --threading=omp --mode=release"
	echo "$SCONS"
	$SCONS || exit
done

