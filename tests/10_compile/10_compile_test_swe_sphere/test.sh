#! /bin/bash

cd "$SWEET_ROOT"


echo
SCONS="scons --program=swe_sphere --gui=disable --sphere-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS  || exit

echo
SCONS="scons --program=swe_sphere --gui=disable --sphere-spectral-space=enable --plane-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS  || exit

