#! /bin/bash

SWEETROOT="/glade/u/home/martins/workspace/sweet/benchmarks_sphere/sph_rexi_linear_paper_gaussian_ts_comparison_earth_scale_cheyenne_performance/script_g9.80616_h10000_f7.292e-05_a6371220_u0_U0_fsph0_tsm_l_rexi_tso0_tsob1_C129600_REXITER_m00004192_h0.15_nrm1_hlf0_bf0_ext00_M0128_MPI_space01_time512/../../../"
cd "$SWEETROOT"

scons  --mode=release --compiler=intel --debug-symbols=enable --simd=enable --mic=disable --fortran-source=enable --program-binary-name= --sweet-mpi=disable --threading=off --rexi-thread-parallel-sum=enable --numa-block-allocator=2 --program=swe_sphere --parareal=none --libpfasst=disable --pfasst-cpp=disable --libfft=enable --libsph=enable --mkl=disable --plane-spectral-space=disable --plane-spectral-dealiasing=disable --sphere-spectral-space=enable --sphere-spectral-dealiasing=enable --libxml=disable --compiler-c-exec= --compiler-cpp-exec= --compiler-fortran-exec= --gui=disable --quadmath=enable

