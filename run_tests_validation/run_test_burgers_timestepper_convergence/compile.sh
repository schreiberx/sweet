#! /bin/bash

cd ../../
#. ./local_software/env_vars.sh
sweetenv

scons --program=burgers --gui=disable --plane-spectral-space=enable --parareal=none --plane-spectral-dealiasing=enable --mode=release --threading=omp --rexi-thread-parallel-sum=disable --compiler=gnu -j 4


