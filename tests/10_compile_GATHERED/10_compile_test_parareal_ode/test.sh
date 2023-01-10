#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"

## --plane-spectral-space=disable->enable
echo
echo "PARAREAL ODE"
SCONS="scons --program=parareal_ode --parareal=serial --gui=disable --plane-spectral-space=enable --mode=debug "
echo "$SCONS"
$SCONS || exit

