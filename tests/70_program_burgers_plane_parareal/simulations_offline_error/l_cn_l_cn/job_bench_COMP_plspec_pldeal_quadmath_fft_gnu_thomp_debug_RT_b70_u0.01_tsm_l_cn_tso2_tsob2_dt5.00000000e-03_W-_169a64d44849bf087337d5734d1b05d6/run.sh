#! /bin/bash

#
# Generating function: jobscript_get_header
# Platform: default_gnu
# Job id: RT_b70_u0.01_tsm_l_cn_tso2_tsob2_dt5.00000000e-03_W-00001_M0032_N-001_plansquick_par_5_ptsm_l_cn_pDt_-1_pStore_0
#


# %SCRIPT_HEADER%


# Loading Job environment variables for currently active platform
cd "/home/jcaldass/Development/IME/sweet"
source ./local_software/env_vars.sh "default_gnu" || exit 1

export MPICC=mpicc
export MPICXX=mpic++
export MPIF90=mpif90


SCONS="scons  --mode=debug --compiler=gnu --sanitize= --debug-symbols=enable --simd=enable --mic=disable --fortran-source=disable --lapack=enable --program-binary-name= --sweet-mpi=disable --threading=omp --rexi-thread-parallel-sum=disable --benchmark-timings=disable --rexi-timings-additional-barriers=disable --rexi-allreduce=disable --program=burgers --parareal=serial --parareal-scalar=disable --parareal-plane=disable --parareal-sphere=disable --parareal-plane-swe=disable --parareal-plane-burgers=disable --libpfasst=disable --eigen=disable --libfft=enable --libsph=disable --mkl=disable --plane-spectral-space=enable --plane-spectral-dealiasing=enable --sphere-spectral-space=disable --sphere-spectral-dealiasing=disable --libxml=disable --gui=disable --quadmath=enable -j 4"
echo "$SCONS"
$SCONS || exit 1


# chdir to execution directory
cd "/home/jcaldass/Development/IME/sweet/tests/70_program_burgers_plane_parareal/job_bench_COMP_plspec_pldeal_quadmath_fft_gnu_thomp_debug_RT_b70_u0.01_tsm_l_cn_tso2_tsob2_dt5.00000000e-03_W-_169a64d44849bf087337d5734d1b05d6"

# %SCRIPT_EXEC_PREFIX%


#
# Generating function: jobscript_get_exec_command
# Platform: default_gnu
# Job id: RT_b70_u0.01_tsm_l_cn_tso2_tsob2_dt5.00000000e-03_W-00001_M0032_N-001_plansquick_par_5_ptsm_l_cn_pDt_-1_pStore_0
#


# mpiexec ... would be here without a line break
EXEC="/home/jcaldass/Development/IME/sweet/build/burgers_COMP_plspec_pldeal_quadmath_fft_gnu_thomp_debug  -G 0 -M 32 -N -1 --space-grid-use-c-staggering=0 -S 1 --benchmark-name=70 -v 2 --dt=0.005 -o 0.1 -u 0.01 -t 1.0 --max-wallclock-time -1 -d 12 --timestepping-method=l_cn --timestepping-order=2 --timestepping-order2=2 --rexi-method=file --rexi-sphere-preallocation=0 --rexi-files= --semi-lagrangian-approximate-sphere-geometry=0 --compute-error=0 --reuse-plans=quick --parareal-enable=1 --parareal-coarse-slices=5 --parareal-convergence-threshold=-1 --parareal-verbosity=6 --parareal-max-simulation-time=1.0 --parareal-coarse-timestepping-method=l_cn --parareal-coarse-timestepping-order=2 --parareal-coarse-timestepping-order2=2 --parareal-coarse-timestep-size=-1 --parareal-load-ref-csv-files=1 --parareal-path-ref-csv-files=../simulations_offline_error/l_cn_l_cn/job_benchref_COMP_plspec_pldeal_quadmath_fft_gnu_thomp_debug_RT_b70_u0.01_tsm_ln_erk_tso2_tsob2_dt1.00000000e-_5f04dce99dfc720226b8553e3873afd1 --parareal-load-fine-csv-files=1 --parareal-path-fine-csv-files=../simulations_offline_error/l_cn_l_cn/job_bench_COMP_plspec_pldeal_quadmath_fft_gnu_thomp_debug_RT_b70_u0.01_tsm_l_cn_tso2_tsob2_dt5.00000000e-03_W-_1f9fffb0d26250a81ef633321945923c --parareal-store-iterations=0"
echo "$EXEC"
$EXEC || exit 1

# %SCRIPT_EXEC_SUFFIX%


#
# Generating function: jobscript_get_footer
# Platform: default_gnu
# Job id: RT_b70_u0.01_tsm_l_cn_tso2_tsob2_dt5.00000000e-03_W-00001_M0032_N-001_plansquick_par_5_ptsm_l_cn_pDt_-1_pStore_0
#


# %SCRIPT_FOOTER%
return 2>/dev/null; exit 0