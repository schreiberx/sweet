#! /bin/bash

SWEETROOT="/home/hpc/pr63so/di69fol/workspace/sweet/benchmarks_plane/polvani/script_swe_plane_polvani_N_bpolvani_g9.81_h10000_f7.2921e-05_p0_a6371220_u2e-05_U8_rob1_fsph0_tsm_ln_erk_tso4_tsob1_C000.0025_PR0.4_PF0.1_M0200_MPI_space01_time001/../../../"
cd "$SWEETROOT"

scons  --mode=release --compiler=gcc --debug-symbols=enable --simd=enable --mic=disable --fortran-source=disable --program-binary-name= --sweet-mpi=disable --threading=omp --rexi-thread-parallel-sum=disable --numa-block-allocator=2 --program=swe_plane --parareal=none --libpfasst=disable --pfasst-cpp=disable --libfft=enable --libsph=disable --mkl=disable --plane-spectral-space=enable --plane-spectral-dealiasing=enable --sphere-spectral-space=disable --sphere-spectral-dealiasing=disable --libxml=disable --compiler-c-exec= --compiler-cpp-exec= --compiler-fortran-exec= --gui=disable --quadmath=enable

