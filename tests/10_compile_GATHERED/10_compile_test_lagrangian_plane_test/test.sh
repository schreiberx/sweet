#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"

echo
echo "LAGRANGIAN PLANE TEST"
echo

SCONS="scons --program=lagrangian_plane_test --gui=disable --plane-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS || exit

###SCONS="scons --program=lagrangian_plane_test --gui=disable --plane-spectral-space=disable --mode=debug" 
###echo "$SCONS"
###$SCONS || exit

mule.benchmark.cleanup_all || exit 1
