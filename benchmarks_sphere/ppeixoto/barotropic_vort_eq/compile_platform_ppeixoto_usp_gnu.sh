#! /bin/bash

# Loading Job environment variables
cd "/home/pedrosp/Reps/sweet"
source ./local_software/env_vars.sh "/home/pedrosp/Reps/sweet/mule/platforms/50_ppeixoto_usp_gnu/env_vars.sh" || exit 1



SCONS="scons  --mode=release --compiler=gnu --debug-symbols=enable --simd=enable --mic=disable --fortran-source=enable --lapack=enable --program-binary-name= --sweet-mpi=enable --threading=omp --rexi-thread-parallel-sum=disable --benchmark-timings=enable --rexi-timings-additional-barriers=disable --rexi-allreduce=disable --numa-block-allocator=2 --program=swe_sphere --parareal=none --libpfasst=disable --eigen=disable --libfft=enable --libsph=enable --mkl=disable --plane-spectral-space=disable --plane-spectral-dealiasing=disable --sphere-spectral-space=enable --sphere-spectral-dealiasing=enable --libxml=disable --gui=disable --quadmath=enable -j 4"
echo "$SCONS"
$SCONS || exit 1
