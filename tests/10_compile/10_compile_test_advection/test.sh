#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"

echo
echo "ADVECTION"


SCONS="scons --program=advection --plane-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS || exit

#### flag SWETT_USE_LIBFFT always required
###SCONS="scons --program=advection --plane-spectral-space=disable --mode=debug"
###echo "$SCONS"
###$SCONS || exit

