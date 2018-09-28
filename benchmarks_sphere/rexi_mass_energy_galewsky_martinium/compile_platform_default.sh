#! /bin/bash

# Loading SWEET environment variables
cd "/home/martin/workspace/sweet"
source ./local_software/env_vars.sh "/home/martin/workspace/sweet/platforms/99_default/env_vars.sh" || exit 1



SCONS="scons  --mode=release --compiler=gnu --debug-symbols=enable --simd=enable --mic=disable --fortran-source=enable --program-binary-name= --sweet-mpi=disable --threading=off --rexi-thread-parallel-sum=enable --rexi-timings=disable --rexi-timings-additional-barriers=disable --numa-block-allocator=2 --program=swe_sphere --parareal=none --libpfasst=disable --pfasst-cpp=disable --eigen=disable --libfft=enable --libsph=enable --mkl=disable --plane-spectral-space=disable --plane-spectral-dealiasing=disable --sphere-spectral-space=enable --sphere-spectral-dealiasing=enable --libxml=disable --compiler-c-exec= --compiler-cpp-exec= --compiler-fortran-exec= --gui=disable --quadmath=disable -j 4"
echo "$SCONS"
$SCONS
