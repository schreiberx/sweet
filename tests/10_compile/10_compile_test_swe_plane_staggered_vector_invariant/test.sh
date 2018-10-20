#! /bin/bash

cd "$SWEET_ROOT"

echo
echo "SWE staggered_vector_invariant"
echo

SCONS="scons --program=swe_plane_staggered_vector_invariant --gui=disable --plane-spectral-space=enable --mode=debug"
echo "$SCONS"
$SCONS  || exit

SCONS="scons --program=swe_plane_staggered_vector_invariant --gui=disable --plane-spectral-space=disable --mode=debug"
echo "$SCONS"
$SCONS  || exit


