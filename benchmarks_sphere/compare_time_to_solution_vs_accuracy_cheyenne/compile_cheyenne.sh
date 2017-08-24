#! /bin/bash

SWEETROOT="/home/martin/workspace/sweet/benchmarks_sphere/compare_time_to_solution_vs_accuracy_cheyenne/script_ln2_bench4_g1_h100000_f7.2921e-05_pde0_a6371220_u0.0_rob1_fsph0_tsm_ln_etdrk_tso2_tsob2_REXICI_n00000128_sx100.0_sy100.0_mu0_primcircle_norm0_half1_preal1_extm00_C002048.0_M0064_MPI_space1_time128/../../../"
cd "$SWEETROOT"

scons  --mode=release --compiler=gnu --debug-symbols=enable --simd=enable --mic=disable --program-binary-name= --sweet-mpi=disable --threading=off --rexi-thread-parallel-sum=enable --numa-block-allocator=2 --program=swe_sphere_rexi --parareal=none --libpfasst=disable --pfasst-cpp=disable --libfft=enable --libsph=enable --mkl=disable --plane-spectral-space=disable --plane-spectral-dealiasing=disable --sphere-spectral-space=enable --sphere-spectral-dealiasing=enable --libxml=disable --gui=disable

