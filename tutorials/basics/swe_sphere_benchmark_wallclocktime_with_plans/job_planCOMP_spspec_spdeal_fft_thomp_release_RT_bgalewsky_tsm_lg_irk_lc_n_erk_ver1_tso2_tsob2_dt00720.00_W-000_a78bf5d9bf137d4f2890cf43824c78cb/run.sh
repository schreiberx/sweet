#! /bin/bash

#
# Generating function: jobscript_get_header
# Platform: default_gcc
# Job id: RT_bgalewsky_tsm_lg_irk_lc_n_erk_ver1_tso2_tsob2_dt00720.00_W-00001_M0064_planssave
#



export OMP_NUM_THREADS=8
export OMP_DISPLAY_ENV=VERBOSE
# %SCRIPT_HEADER%


# Loading Job environment variables for currently active platform
cd "/home/martin/workspace/sweet_REPOS"
source ./local_software/env_vars.sh "default_gcc" || exit 1

export MPICC=mpicc
export MPICXX=mpic++
export MPIF90=mpif90


SCONS="scons  --mode=release --debug-symbols=enable --simd=enable --fortran-source=enable --lapack=enable --program-binary-name= --sweet-mpi=disable --threading=omp --rexi-thread-parallel-sum=disable --benchmark-timings=disable --rexi-timings-additional-barriers=disable --rexi-allreduce=disable --program=swe_sphere --parareal=none --parareal-scalar=disable --parareal-plane=disable --parareal-sphere=disable --parareal-plane-swe=disable --parareal-plane-burgers=disable --xbraid=none --xbraid-scalar=disable --xbraid-plane=disable --xbraid-sphere=disable --xbraid-plane-swe=disable --xbraid-plane-burgers=disable --libpfasst=disable --eigen=disable --libfft=enable --libsph=enable --mkl=disable --plane-spectral-space=disable --plane-spectral-dealiasing=disable --sphere-spectral-space=enable --sphere-spectral-dealiasing=enable --libxml=disable --gui=disable --quadmath=disable -j 4"
echo "$SCONS"
$SCONS || exit 1


# chdir to execution directory
cd "/home/martin/workspace/sweet_REPOS/tutorials/basics/swe_sphere_benchmark_wallclocktime_with_plans/job_planCOMP_spspec_spdeal_fft_thomp_release_RT_bgalewsky_tsm_lg_irk_lc_n_erk_ver1_tso2_tsob2_dt00720.00_W-000_a78bf5d9bf137d4f2890cf43824c78cb"

# %SCRIPT_EXEC_PREFIX%


#
# Generating function: jobscript_get_exec_command
# Platform: default_gcc
# Job id: RT_bgalewsky_tsm_lg_irk_lc_n_erk_ver1_tso2_tsob2_dt00720.00_W-00001_M0064_planssave
#


# mpiexec ... would be here without a line break
EXEC="/home/martin/workspace/sweet_REPOS/build/swe_sphere_COMP_spspec_spdeal_fft_thomp_release  -M 64 --space-grid-use-c-staggering=0 -S 1 --benchmark-name=galewsky -v 0 --dt=720 -o -1 --output-file-name=- -t 0 --max-wallclock-time -1 --instability-checks=1 -d 12 --timestepping-method=lg_irk_lc_n_erk_ver1 --timestepping-order=2 --timestepping-order2=2 --rexi-method=file --rexi-sphere-preallocation=0 --rexi-files= --semi-lagrangian-approximate-sphere-geometry=0 --compute-error=0 --num-threads-space=-1 --reuse-plans=save"
echo "$EXEC"
$EXEC || exit 1

# %SCRIPT_EXEC_SUFFIX%

cp ./shtns_cfg ../ 2>/dev/null
cp ./shtns_cfg_fftw ../ 2>/dev/null



#
# Generating function: jobscript_get_footer
# Platform: default_gcc
# Job id: RT_bgalewsky_tsm_lg_irk_lc_n_erk_ver1_tso2_tsob2_dt00720.00_W-00001_M0064_planssave
#


# %SCRIPT_FOOTER%
return 2>/dev/null; exit 0