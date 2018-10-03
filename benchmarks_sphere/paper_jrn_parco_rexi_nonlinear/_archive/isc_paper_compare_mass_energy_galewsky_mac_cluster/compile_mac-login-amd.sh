#! /bin/bash

SWEETROOT="/home/hpc/pr63so/di69fol/workspace/sweet/benchmarks_sphere/isc_paper_compare_mass_energy_galewsky_mac_cluster/script_ln2_g9.81_h10000_f7.2921e-05_a6371220_u0.0_U0_fsph0_tsm_l_rexi_n_etdrk_tso2_tsob2_C000720_REXICI_n00000128_mr10.0_mi30.0_prcircle_gfs0.0000E+00_gfd0.0000E+00_gfe0.0000E+00_nrm0_hlf0_bf0_ext00_M0128_MPI_space01_time001/../../../"
cd "$SWEETROOT"

scons  --mode=release --compiler=gnu --debug-symbols=enable --simd=enable --mic=disable --fortran-source=enable --program-binary-name= --sweet-mpi=disable --threading=off --rexi-thread-parallel-sum=enable --rexi-timings=disable --numa-block-allocator=2 --program=swe_sphere --parareal=none --libpfasst=disable --pfasst-cpp=disable --eigen=disable --libfft=enable --libsph=enable --mkl=disable --plane-spectral-space=disable --plane-spectral-dealiasing=disable --sphere-spectral-space=enable --sphere-spectral-dealiasing=enable --libxml=disable --compiler-c-exec= --compiler-cpp-exec= --compiler-fortran-exec= --gui=disable --quadmath=enable

