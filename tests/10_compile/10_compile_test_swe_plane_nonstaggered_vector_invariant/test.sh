#! /bin/bash

cd "$SWEET_ROOT"

echo
echo "SWE nonstaggered_vector_invariant"
echo

SCONS="scons --program=swe_plane_nonstaggered_vector_invariant --gui=disable --mode=debug"
echo "$SCONS"
$SCONS  || exit

