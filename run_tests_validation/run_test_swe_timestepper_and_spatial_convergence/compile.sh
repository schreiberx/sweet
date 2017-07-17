#! /bin/bash

cd ../../
. ./local_software/env_vars.sh

scons --program=swe_plane_rexi --gui=disable --plane-spectral-space=enable --mode=release --threading=omp --rexi-thread-parallel-sum=disable --compiler=gnu -j 4

scons --program=swe_plane_rexi --gui=disable --plane-spectral-space=disable --mode=release --threading=omp --rexi-thread-parallel-sum=disable --compiler=gnu -j 4

