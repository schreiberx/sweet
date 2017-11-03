#! /bin/bash

cd ../../
. ./local_software/env_vars.sh

scons --program=swe_plane --gui=disable --plane-spectral-space=enable --mode=release --threading=omp --rexi-thread-parallel-sum=disable --compiler=gnu -j 4


