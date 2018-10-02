#! /bin/bash

cd "$SWEET_ROOT"

for i in 0 1 2 3; do
	SCONS="scons --program=advection --numa-block-allocator=$i --threading=omp --mode=debug"
	echo "$SCONS"
	$SCONS || exit
done

echo
echo "ADVECTION"

SCONS="scons --program=advection --plane-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS || exit

SCONS="scons --program=advection --plane-spectral-space=disable --mode=debug"
echo "$SCONS"
$SCONS || exit

