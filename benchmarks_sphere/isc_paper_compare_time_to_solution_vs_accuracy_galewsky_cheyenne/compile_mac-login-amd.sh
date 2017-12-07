#! /bin/bash

SWEETROOT="/home/hpc/pr63so/di69fol/workspace/sweet/benchmarks_sphere/isc_paper_compare_time_to_solution_vs_accuracy_galewsky_mac_cluster/script_ln2_b100_g9.81_h10000_f7.2921e-05_a6371220_u0.0_U0_fsph0_tsm_l_rexi_n_erk_ver1_tso2_tsob2_C0920_REXICI_n00000256_mr5.0_mi60.0_prcircle_gfs5.0000E+01_gfd1.3000E+02_gfe1.0000E+01_nrm0_hlf0_bf1e-16_ext00_M0128_MPI_space01_time256/../../../"
cd "$SWEETROOT"

scons  --mode=release --compiler=gnu --debug-symbols=enable --simd=enable --mic=disable --fortran-source=enable --program-binary-name= --sweet-mpi=disable --threading=off --rexi-thread-parallel-sum=enable --numa-block-allocator=2 --program=swe_sphere --parareal=none --libpfasst=disable --pfasst-cpp=disable --libfft=enable --libsph=enable --mkl=disable --plane-spectral-space=disable --plane-spectral-dealiasing=disable --sphere-spectral-space=enable --sphere-spectral-dealiasing=enable --libxml=disable --compiler-c-exec= --compiler-cpp-exec= --compiler-fortran-exec= --gui=disable --quadmath=enable

