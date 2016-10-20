#! /bin/bash

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

echo
echo "BURGERS"

SCONS="scons --program=burgers --parareal=serial --plane-spectral-space=enable --mode=debug"
echo $SCONS
$SCONS || exit

echo
echo "LAGRANGIAN TEST"
SCONS="scons --program=lagrangian_test --gui=disable --plane-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS || exit
SCONS="scons --program=lagrangian_test --gui=disable --plane-spectral-space=disable --mode=debug" 
echo "$SCONS"
$SCONS || exit

echo
echo "PARAREAL ODE"
SCONS="scons --program=parareal_ode --parareal=serial --gui=disable --plane-spectral-space=disable --mode=debug "
echo "$SCONS"
$SCONS || exit

echo
echo "SPECTRAL VISUALIZATION"
SCONS="scons --program=spectral_visualization --gui=enable --plane-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS  || exit


echo
echo "SPHERICAL HARMONICS OUTPUT"
SCONS="scons --program=sph_output_spherical_harmonics --gui=disable --sphere-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS  || exit


echo
echo "SWE nonstaggered_advective"
SCONS="scons --program=swe_nonstaggered_advective --gui=enable --mode=debug"
echo "$SCONS"
$SCONS  || exit
SCONS="scons --program=swe_nonstaggered_advective --plane-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS  || exit


echo
echo "SWE nonstaggered_advective_linear_only"
SCONS="scons --program=swe_nonstaggered_advective_linear_only --plane-spectral-space=disable --mode=debug"
echo "$SCONS"
$SCONS  || exit
SCONS="scons --program=swe_nonstaggered_advective_linear_only --plane-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS  || exit


echo
echo "SWE nonstaggered_vector_invariant"
SCONS="scons --program=swe_nonstaggered_vector_invariant --gui=enable --mode=debug"
echo "$SCONS"
$SCONS  || exit

echo
echo "SWE REXI"
SCONS="scons --program=swe_rexi --gui=disable --plane-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS  || exit
SCONS="scons --program=swe_rexi --gui=disable --plane-spectral-space=disable --mode=debug"
echo "$SCONS"
$SCONS  || exit
SCONS="scons --program=swe_rexi --gui=disable --plane-spectral-space=enable --plane-spectral-dealiasing=disable --mode=debug"
echo "$SCONS"
$SCONS  || exit

echo
echo "SPHERICAL HARMONICS SWE AND REXI"
SCONS="scons --program=swe_sph_and_rexi --gui=disable --plane-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS  || exit

echo
echo "SWE staggered_vector_invariant"
SCONS="scons --program=swe_staggered_vector_invariant --gui=enable --plane-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS  || exit
SCONS="scons --program=swe_staggered_vector_invariant --gui=enable --plane-spectral-space=disable --mode=debug"
echo "$SCONS"
$SCONS  || exit

