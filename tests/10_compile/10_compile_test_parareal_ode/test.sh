#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"

## --plane-spectral-space=disable->enable
echo
echo "PARAREAL ODE"
##SCONS="scons --program=parareal_ode --parareal=serial --gui=disable --plane-spectral-space=enable --mode=debug --scalar-type=complex --N-ode=1"
SCONS="scons --program=parareal_ode --parareal=none --gui=disable --plane-spectral-space=enable --mode=debug --scalar-type=complex --N-ode=1 --eigen=enable"
echo "$SCONS"
$SCONS || exit

