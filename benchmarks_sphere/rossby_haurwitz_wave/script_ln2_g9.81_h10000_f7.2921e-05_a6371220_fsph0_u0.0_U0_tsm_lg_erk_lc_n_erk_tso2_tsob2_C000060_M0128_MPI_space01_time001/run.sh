#!/bin/bash




cd "/home/martin/workspace/sweet/benchmarks_sphere/rossby_haurwitz_wave/script_ln2_g9.81_h10000_f7.2921e-05_a6371220_fsph0_u0.0_U0_tsm_lg_erk_lc_n_erk_tso2_tsob2_C000060_M0128_MPI_space01_time001"

BASEDIR="`pwd`"
#rm -f ./prog_h_*
#rm -f ./prog_u_*
#rm -f ./prog_v_*

SWEETROOT="/home/martin/workspace/sweet/benchmarks_sphere/rossby_haurwitz_wave/script_ln2_g9.81_h10000_f7.2921e-05_a6371220_fsph0_u0.0_U0_tsm_lg_erk_lc_n_erk_tso2_tsob2_C000060_M0128_MPI_space01_time001/../../../"
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


scons  --mode=release --compiler=gnu --debug-symbols=enable --simd=enable --mic=disable --fortran-source=enable --program-binary-name= --sweet-mpi=disable --threading=off --rexi-thread-parallel-sum=disable --rexi-timings=disable --rexi-timings-additional-barriers=disable --numa-block-allocator=2 --program=swe_sphere --parareal=none --libpfasst=disable --pfasst-cpp=disable --eigen=disable --libfft=enable --libsph=enable --mkl=disable --plane-spectral-space=disable --plane-spectral-dealiasing=disable --sphere-spectral-space=enable --sphere-spectral-dealiasing=enable --libxml=disable --compiler-c-exec= --compiler-cpp-exec= --compiler-fortran-exec= --gui=disable --quadmath=enable



cd "$BASEDIR"
pwd

EXEC="$SWEETROOT/build/swe_sphere_spspec_spdeal_quadmath_numa2_fft_gnu_release  -g 9.81 -H 10000 -f 7.2921e-05 -F 0 -a 6371220 -M 128 --pde-id 0 --staggering=0 -S 1 -X 1 -Y 1 -s 4 --benchmark=rossby_haurwitz_wave -v 2 --dt=60 -o 3600 -u 0.0 -t 3600 --stability-checks=0 --use-linear-div=0 --timestepping-method=lg_erk_lc_n_erk --timestepping-order=2 --timestepping-order2=2 --normal-mode-analysis-generation=0 --rexi-method= --polvani-rossby=-1.0 --polvani-froude=-1.0 --use-robert-functions=1 --compute-error=0"


echo "$EXEC"
pwd
#ln -s "$SWEETROOT/data/" "$BASEDIR/data"   #Symlink for GUI directory, if necessary
$EXEC || exit 1
