#! /bin/bash

cd "$SWEET_ROOT"

echo
echo "ADVECTION"


SCONS="scons --program=advection --plane-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS || exit

SCONS="scons --program=advection --plane-spectral-space=disable --mode=debug"
echo "$SCONS"
$SCONS || exit

