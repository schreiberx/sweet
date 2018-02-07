#! /bin/bash

SWEETROOT="/home/hpc/pr63so/uq_weather_nonrode/sweet/benchmarks_sphere/uq_galewsky_mac_cluster/script_ln2_bgu8.8000E+01_bgh1.3200E+02_bgp8.6394E-01_g10.786776_h11000.0_f8.0212e-05_a6371220_u0.0_fsph0_tsm_lg_irk_lc_n_erk_ver1_tso2_tsob2_C000010_M0128_MPI_space01_time001/../../../"
cd "$SWEETROOT"

scons  --mode=release --compiler=gnu --debug-symbols=enable --simd=enable --mic=disable --fortran-source=enable --program-binary-name= --sweet-mpi=disable --threading=off --rexi-thread-parallel-sum=enable --rexi-timings=disable --rexi-timings-additional-barriers=disable --numa-block-allocator=2 --program=swe_sphere --parareal=none --libpfasst=disable --pfasst-cpp=disable --eigen=disable --libfft=enable --libsph=enable --mkl=disable --plane-spectral-space=disable --plane-spectral-dealiasing=disable --sphere-spectral-space=enable --sphere-spectral-dealiasing=enable --libxml=disable --compiler-c-exec= --compiler-cpp-exec= --compiler-fortran-exec= --gui=disable --quadmath=enable

