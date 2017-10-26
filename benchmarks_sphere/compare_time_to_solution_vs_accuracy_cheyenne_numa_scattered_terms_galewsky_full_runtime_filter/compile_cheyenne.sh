#! /bin/bash

SWEETROOT="/home/martin/workspace/sweet/benchmarks_sphere/compare_time_to_solution_vs_accuracy_cheyenne_numa_scattered_terms_galewsky_full_runtime_filter/script_ln2_b100_g9.81_h10000_f7.2921e-05_p0_a6371220_u0.0_U0_rob1_fsph0_tsm_lg_rexi_lc_n_erk_tso2_tsob2_REXICI_n00000128_mr10.0_mi30.0_prcircle_gf0.0001_nrm0_hlf0_pre1_ext00_C00002560_M0128_MPI_space01_time128/../../../"
cd "$SWEETROOT"

scons  --mode=release --compiler=gnu --debug-symbols=enable --simd=enable --mic=disable --fortran-source=enable --program-binary-name= --sweet-mpi=disable --threading=off --rexi-thread-parallel-sum=enable --numa-block-allocator=2 --program=swe_sphere_rexi --parareal=none --libpfasst=disable --pfasst-cpp=disable --libfft=enable --libsph=enable --mkl=disable --plane-spectral-space=disable --plane-spectral-dealiasing=disable --sphere-spectral-space=enable --sphere-spectral-dealiasing=enable --libxml=disable --compiler-c-exec= --compiler-cpp-exec= --compiler-fortran-exec= --gui=disable --quadmath=enable

