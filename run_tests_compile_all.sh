#! /bin/bash

source ./local_software/env_vars.sh

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


if [ "x" != "x$DISPLAY" ]; then 
	echo
	echo "SPECTRAL VISUALIZATION"
	SCONS="scons --program=spectral_visualization --gui=enable --plane-spectral-space=enable --mode=debug"
	echo "$SCONS"
	$SCONS  || exit
fi


echo
echo "SPHERICAL HARMONICS OUTPUT"
SCONS="scons --program=sphere_output_spherical_harmonics --gui=disable --sphere-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS  || exit



echo
echo "SWE nonstaggered_vector_invariant"
SCONS="scons --program=swe_plane_nonstaggered_vector_invariant --gui=disable --mode=debug"
echo "$SCONS"
$SCONS  || exit

echo
echo "SWE REXI"
SCONS="scons --program=swe_plane --gui=disable --plane-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS  || exit

echo
echo "SWE REXI"
SCONS="scons --program=swe_plane --gui=disable --plane-spectral-space=disable --libfft=enable --mode=debug"
echo "$SCONS"
$SCONS  || exit

echo
echo "SWE REXI"
SCONS="scons --program=swe_plane --gui=disable --plane-spectral-space=disable --libfft=enable --mode=debug"
echo "$SCONS"
$SCONS  || exit


SCONS="scons --program=swe_plane --gui=disable --plane-spectral-space=enable --plane-spectral-dealiasing=disable --mode=debug"
echo "$SCONS"
$SCONS  || exit



mpiCC -v 2>&1 2> /dev/null
if [ $? -eq 0 ]; then
	SCONS="scons --program=swe_plane --sweet-mpi=enable --rexi-thread-parallel-sum=enable --threading=off"
	echo "$SCONS"
	$SCONS  || exit
fi



echo
echo "SWE staggered_vector_invariant"
SCONS="scons --program=swe_plane_staggered_vector_invariant --gui=disable --plane-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS  || exit
SCONS="scons --program=swe_plane_staggered_vector_invariant --gui=disable --plane-spectral-space=disable --mode=debug"
echo "$SCONS"
$SCONS  || exit


echo
echo "SPHERICAL HARMONICS SWE AND REXI"
SCONS="scons --program=swe_sphere --gui=disable --sphere-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS  || exit


echo "********************************************"
echo "******************** FIN *******************"
echo "********************************************"
