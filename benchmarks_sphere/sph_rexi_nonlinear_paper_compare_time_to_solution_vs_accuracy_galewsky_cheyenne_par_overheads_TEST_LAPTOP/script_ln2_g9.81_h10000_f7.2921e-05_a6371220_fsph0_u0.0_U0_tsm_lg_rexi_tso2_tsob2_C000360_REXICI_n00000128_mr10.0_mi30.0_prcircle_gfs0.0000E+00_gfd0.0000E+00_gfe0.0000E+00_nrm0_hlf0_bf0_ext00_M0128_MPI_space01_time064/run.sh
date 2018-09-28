#!/bin/bash




cd "/home/martin/workspace/sweet/benchmarks_sphere/sph_rexi_nonlinear_paper_compare_time_to_solution_vs_accuracy_galewsky_cheyenne_par_overheads_TEST_LAPTOP/script_ln2_g9.81_h10000_f7.2921e-05_a6371220_fsph0_u0.0_U0_tsm_lg_rexi_tso2_tsob2_C000360_REXICI_n00000128_mr10.0_mi30.0_prcircle_gfs0.0000E+00_gfd0.0000E+00_gfe0.0000E+00_nrm0_hlf0_bf0_ext00_M0128_MPI_space01_time064"

BASEDIR="`pwd`"
#rm -f ./prog_h_*
#rm -f ./prog_u_*
#rm -f ./prog_v_*

SWEETROOT="/home/martin/workspace/sweet/benchmarks_sphere/sph_rexi_nonlinear_paper_compare_time_to_solution_vs_accuracy_galewsky_cheyenne_par_overheads_TEST_LAPTOP/script_ln2_g9.81_h10000_f7.2921e-05_a6371220_fsph0_u0.0_U0_tsm_lg_rexi_tso2_tsob2_C000360_REXICI_n00000128_mr10.0_mi30.0_prcircle_gfs0.0000E+00_gfd0.0000E+00_gfe0.0000E+00_nrm0_hlf0_bf0_ext00_M0128_MPI_space01_time064/../../../"
cd "$SWEETROOT"

pwd

# Always load local software
#is this really the root?
if test -e ./local_software/env_vars.sh ; then
	source ./local_software/env_vars.sh || exit 1
else #try ../
	echo "Warning: changing SWEETROOT directory"	
	cd ..
	SWEETROOT="`pwd`"
	pwd
	source ./local_software/env_vars.sh || exit 1
fi
#source ./local_software/env_vars.sh || exit 1

#make clean || exit 1


SCONS="scons  --mode=release --compiler=gnu --debug-symbols=enable --simd=enable --mic=disable --fortran-source=enable --program-binary-name= --sweet-mpi=enable --threading=off --rexi-thread-parallel-sum=disable --rexi-timings=enable --rexi-timings-additional-barriers=enable --numa-block-allocator=2 --program=swe_sphere --parareal=none --libpfasst=disable --pfasst-cpp=disable --eigen=disable --libfft=enable --libsph=enable --mkl=disable --plane-spectral-space=disable --plane-spectral-dealiasing=disable --sphere-spectral-space=enable --sphere-spectral-dealiasing=enable --libxml=disable --compiler-c-exec= --compiler-cpp-exec= --compiler-fortran-exec= --gui=disable --quadmath=disable -j 4"
echo "$SCONS"
$SCONS || exit 1



cd "$BASEDIR"
pwd

EXEC="$SWEETROOT/build/swe_sphere_spspec_spdeal_mpi_rxtime_rxtbar_numa2_fft_gnu_release  -g 9.81 -H 10000 -f 7.2921e-05 -F 0 -a 6371220 -M 128 --pde-id 0 --staggering=0 -S 1 -X 1 -Y 1 -s 100 --benchmark= -v 0 --dt=360 -o 432000 -O - -u 0.0 -t 432000 --stability-checks=1 --use-linear-div=0 --use-local-visc=0 --timestepping-method=lg_rexi --timestepping-order=2 --timestepping-order2=2 --normal-mode-analysis-generation=0 --rexi-method=ci --rexi-half=0 --rexi-normalization=0 --rexi-sphere-preallocation=1 --rexi-use-direct-solution=0 --rexi-ext-modes=0 --rexi-ci-n=128 --rexi-ci-max-real=10.0 --rexi-ci-max-imag=30.0 --rexi-ci-primitive=circle --rexi-ci-gaussian-filter-scale=0.0 --rexi-ci-gaussian-filter-dt-norm=0.0 --rexi-ci-gaussian-filter-exp-N=0.0 --polvani-rossby=-1.0 --polvani-froude=-1.0 --use-robert-functions=1 --compute-error=0 --shtns-use-plans=0"


echo "$EXEC"
pwd
#ln -s "$SWEETROOT/data/" "$BASEDIR/data"   #Symlink for GUI directory, if necessary
$EXEC || exit 1
