#! /bin/bash

SWEETROOT="/glade/u/home/martins/workspace/sweet/benchmarks_sphere/compare_time_to_solution_vs_accuracy_cheyenne/script_ln2_b4_g1_h100000_f7.2921e-05_p0_a6371220_u0.0_rob1_fsph0_tsm_ln_etdrk_tso2_tsob2_REXICI_n00000064_sx50.0_sy50.0_mu0_prcircle_nrm0_hlf1_pre1_ext00_C002048.0_M0128_MPI_space1_time64/../../../"
cd "$SWEETROOT"

scons  --mode=release --compiler=intel --debug-symbols=enable --simd=enable --mic=disable --fortran-source=enable --program-binary-name= --sweet-mpi=enable --threading=off --rexi-thread-parallel-sum=enable --numa-block-allocator=0 --program=swe_sphere_rexi --parareal=none --libpfasst=disable --pfasst-cpp=disable --libfft=enable --libsph=enable --mkl=disable --plane-spectral-space=disable --plane-spectral-dealiasing=disable --sphere-spectral-space=enable --sphere-spectral-dealiasing=enable --libxml=disable --compiler-c-exec=mpicc --compiler-cpp-exec=mpicxx --compiler-fortran-exec=mpif90 --gui=disable

