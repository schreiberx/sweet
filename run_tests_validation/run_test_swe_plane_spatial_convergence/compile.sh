#! /bin/bash

cd ../../
. ./local_software/env_vars.sh

#spectral methods without rexi
scons --program=swe_plane_rexi --gui=disable --plane-spectral-space=enable --libfft=enable --mode=release --threading=omp --rexi-thread-parallel-sum=disable --compiler=gnu -j 4

#Finite differences mode
scons --program=swe_plane_rexi --gui=disable --plane-spectral-space=disable --libfft=enable --mode=release --threading=omp --rexi-thread-parallel-sum=disable --compiler=gnu --plane-spectral-dealiasing=disable -j 4 


scons --program=swe_plane_rexi --gui=disable --plane-spectral-space=enable --libfft=enable --mode=release --threading=omp --rexi-thread-parallel-sum=disable --plane-spectral-dealiasing=disable --compiler=gnu -j 4
