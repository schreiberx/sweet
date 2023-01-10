#! /bin/bash

cd "$MULE_SOFTWARE_ROOT"

echo
echo "SPHERICAL HARMONICS OUTPUT"
echo

SCONS="scons --program=sphere_output_spherical_harmonics --gui=disable --sphere-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS  || exit


