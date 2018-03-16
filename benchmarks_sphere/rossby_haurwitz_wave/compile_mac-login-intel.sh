#! /bin/bash

SWEETROOT="/home/hpc/pr63so/di69fol/workspace/sweet_intel/benchmarks_sphere/rossby_haurwitz_wave/script_ln2_g9.81_h10000_f7.2921e-05_a6371220_fsph0_u0.0_U0_tsm_lg_erk_lc_n_erk_tso2_tsob2_C000512_M0128_MPI_space16_time001/../../../"
cd "$SWEETROOT"

scons  --mode=release --compiler=gnu --debug-symbols=enable --simd=enable --mic=disable --fortran-source=enable --program-binary-name= --sweet-mpi=disable --threading=off --rexi-thread-parallel-sum=disable --rexi-timings=disable --rexi-timings-additional-barriers=disable --numa-block-allocator=2 --program=swe_sphere --parareal=none --libpfasst=disable --pfasst-cpp=disable --eigen=disable --libfft=enable --libsph=enable --mkl=disable --plane-spectral-space=disable --plane-spectral-dealiasing=disable --sphere-spectral-space=enable --sphere-spectral-dealiasing=enable --libxml=disable --compiler-c-exec= --compiler-cpp-exec= --compiler-fortran-exec= --gui=disable --quadmath=enable

