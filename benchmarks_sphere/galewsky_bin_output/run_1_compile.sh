#! /bin/bash

BASEDIR="`pwd`"
SWEETROOT="../../"
cd "$SWEETROOT"

make clean
scons --program=swe_sphere --gui=disable --plane-spectral-space=disable --sphere-spectral-space=enable --threading=omp --mode=release || exit 1

